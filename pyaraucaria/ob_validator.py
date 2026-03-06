from jsonschema import Draft7Validator
from copy import deepcopy
import re
import os
import yaml


class ObsValidator:

    def __init__(self, base_schema: dict, command_rules: dict):
        self.base_schema = base_schema
        self.command_rules = command_rules

    # converts to standard flat dictionary - I need this for purposes
    @staticmethod
    def convert_to_obdict(ob: dict) -> dict:
        result = {"name": None, "ra": None, "dec": None, "command_name": None}

        subcommands = ob.get("subcommands", [])
        if isinstance(subcommands, dict):
            subcommands = [subcommands]

        for sub in subcommands:
            if "command_name" in sub:
                result["command_name"] = sub["command_name"]
            if "kwargs" in sub and isinstance(sub["kwargs"], dict):
                result.update(sub["kwargs"])
            if "args" in sub and isinstance(sub["args"], list):
                if len(sub["args"]) == 1:
                    result["name"] = sub["args"][0]
                elif len(sub["args"]) == 2:
                    result["ra"] = sub["args"][0]
                    result["dec"] = sub["args"][1]
                elif len(sub["args"]) == 3:
                    result["name"] = sub["args"][0]
                    result["ra"] = sub["args"][1]
                    result["dec"] = sub["args"][2]

        return result


    @staticmethod
    def clean_none(obs: dict) -> dict:
        return {k: v for k, v in obs.items() if v is not None}

    # loading schema as dictionary
    @staticmethod
    def load_schema(name: str) -> dict:
        lib_dir = os.path.dirname(__file__)
        schemas_dir = os.path.join(lib_dir, "..", "schemas")
        if ".yaml" not in name:
            name += ".yaml"
        file_path = os.path.join(schemas_dir, name)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Schema not found: {file_path}")
        with open(file_path, "r") as f:
            return yaml.safe_load(f)

    # converts types to declared in schema
    @staticmethod
    def convert_types(obs: dict, schema: dict) -> dict:
        properties = schema.get("properties", {})
        converted = {}

        for key, value in obs.items():
            if key not in properties:
                converted[key] = value
                continue

            prop = properties[key]
            expected = prop.get("type")
            if not expected:
                converted[key] = value
                continue
            if not isinstance(expected, list):
                expected = [expected]

            new_value = value
            converted_successfully = False

            for typ in expected:
                try:
                    if typ == "integer":
                        new_value = int(value)
                        converted_successfully = True
                        break
                    elif typ == "number":
                        new_value = float(value)
                        converted_successfully = True
                        break
                    elif typ == "boolean":
                        if isinstance(value, str):
                            if value.lower() in ["true", "1"]:
                                new_value = True
                                converted_successfully = True
                                break
                            if value.lower() in ["false", "0"]:
                                new_value = False
                                converted_successfully = True
                                break
                        elif isinstance(value, bool):
                            new_value = value
                            converted_successfully = True
                            break
                    elif typ == "string":
                        new_value = str(value)
                        converted_successfully = True
                        break
                    elif typ == "null":
                        if value is None:
                            new_value = None
                            converted_successfully = True
                            break
                except (ValueError, TypeError):
                    pass

            converted[key] = new_value if converted_successfully else value

        return converted

    # walidacja rules
    @staticmethod
    def validate_rules(obs: dict, rules: dict) -> dict:
        result = {}

        # required
        for key in rules.get("required", []):
            if key not in obs:
                result[key] = None

        # one_of - dziwne, ale dziala
        for subschema in rules.get("one_of", []):
            if not any(k not in obs for k in subschema):
                continue
            satisfied = any(all(k in obs for k in s) for s in rules.get("one_of", []))
            if not satisfied:
                for s in rules.get("one_of", []):
                    for k in s:
                        result[k] = False

        # one_of_group
        for group in rules.get("one_of_group", []):
            if not any(k in obs for k in group):
                for k in group:
                    result[k] = False

        return result

    # walidacja seq
    @staticmethod
    def validate_seq(obs: dict, allowed_filters=None) -> dict:
        result = {}
        allowed_filters = allowed_filters or []
        seq_str = obs.get("seq")
        cmd = obs.get("command_name")
        if seq_str is None:
            return result

        def _validate(seq_str):
            seq_str = seq_str.strip()
            match = re.match(r"^(\d+)x\((.+)\)$", seq_str)
            if match:
                repeat_count = int(match.group(1))
                inner = match.group(2)
                if repeat_count < 1:
                    return False
                return _validate(inner)

            exposures = [s.strip() for s in seq_str.split(",")]
            for exp in exposures:
                parts = exp.split("/")
                if len(parts) != 3:
                    return False
                try:
                    n_repeat = int(parts[0])
                    filt = parts[1]
                    exp_time_str = parts[2].replace(",", ".")

                    # dla FOCUS mozliwe exp_time = "a"
                    if cmd == "FOCUS" and exp_time_str.lower() == "a":
                        exp_time = 0  #  placeholder
                    else:
                        exp_time = float(exp_time_str)
                        if exp_time < 0:
                            return False

                    if n_repeat < 1:
                        return False
                    if allowed_filters and filt not in allowed_filters:
                        return False
                except ValueError:
                    return False
            return True

        result["seq"] = _validate(seq_str)
        return result

    def validate_ob(self, obs: dict, overrides: dict = None, allowed_filters=None) -> dict:

        # czysczenie None
        obs_clean = self.clean_none(obs)
        schema = deepcopy(self.base_schema)

        # dopisanie i nadpisanie warunkow intrumentalnych
        if overrides:
            for k, v in overrides.items():
                properties = schema.setdefault("properties", {})
                if k in properties:
                    properties[k].update(v)
                else:
                    properties[k] = v

        properties = schema.get("properties", {})

        # proba konwersji typow
        obs_clean = self.convert_types(obs_clean, schema)

        # walidacja typow
        validator = Draft7Validator(schema)
        result = {key: True for key in obs_clean}
        for error in validator.iter_errors(obs_clean):
            if error.validator in ("required", "oneOf"):
                continue
            if error.path:
                key = error.path[0]
                if key in result:
                    result[key] = False

        # walidacja regul
        cmd_type = obs_clean.get("command_name")
        if cmd_type in self.command_rules:
            result.update(self.validate_rules(obs_clean, self.command_rules[cmd_type]))
        else:
            result["command_name"] = False

        # walidacja seq
        result.update(self.validate_seq(obs_clean, allowed_filters=allowed_filters))

        # dodatkowe zwracanie dla umilenia pracy
        required_fields = []
        allowed_fields = list(properties.keys())
        if cmd_type in self.command_rules:
            rules = self.command_rules[cmd_type]
            required_fields = rules.get("required", [])
            allowed_fields = list(set(list(properties.keys()) + rules.get("allowed", [])))

        valid = all(v is True for v in result.values())
        return {
            "valid": valid,
            "result": result,
            "data": obs_clean,
            "required": required_fields,
            "allowed": allowed_fields
        }










