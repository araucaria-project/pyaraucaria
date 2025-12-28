# **uv Cheat Sheet** tailored for **Poetry** users

### ⚡ Core Philosophy Shift

**Poetry** uses a proprietary config format (`[tool.poetry]`).
**uv** uses standard Python standards (**PEP 621** `[project]`) and explicitly manages the Virtual Environment (`.venv`) in the project root by default.

---

### 🚀 Command Mapping

| Action | Poetry Command | **uv Command** |
| --- | --- | --- |
| **Initialize** | `poetry init` | `uv init` |
| **Add Pkg** | `poetry add numpy` | `uv add numpy` |
| **Dev Pkg** | `poetry add --group dev pytest` | `uv add --dev pytest` |
| **Add Optional** | `poetry add numpy -E plot` | `uv add numpy --optional plot` |
| **Install** | `poetry install` | `uv sync` ✨ *(See Note 1)* |
| **Install Extras** | `poetry install --all-extras` | `uv sync --all-extras` ⚠️ *(See Note 2)* |
| **Remove** | `poetry remove numpy` | `uv remove numpy` |
| **Update** | `poetry update` | `uv lock --upgrade` (lock only)<br>`uv sync --upgrade` (lock + install) |
| **Run** | `poetry run python app.py` | `uv run python app.py` |
| **Shell** | `poetry shell` | *No direct equivalent.* Use `source .venv/bin/activate` |
| **List** | `poetry show` | `uv tree` (visual) or `uv pip list` |

---

### 🐍 Virtual Environments

**Creation**

* **Poetry:** Implicitly creates envs (often centrally located), or you configure `virtualenvs.in-project`.
* **uv:** `uv venv`
* *Note:* uv creates a standard `.venv` folder in your project root by default.



**Activation**

* **Poetry:** `poetry shell`
* **uv:** Standard activation:
* Mac/Linux: `source .venv/bin/activate`
* Windows: `.venv\Scripts\activate`



---

### 📄 pyproject.toml Differences

You will stop using `[tool.poetry]` and start using the standard `[project]` table.

**Poetry Style:**

```toml
[tool.poetry]
name = "my-app"
version = "0.1.0"

[tool.poetry.dependencies]
python = "^3.12"
requests = "^2.31.0"

[tool.poetry.extras]
# Proprietary extras definition
gui = ["PyQt6"]

```

**uv (Standard PEP 621) Style:**

```toml
[project]
name = "my-app"
version = "0.1.0"
requires-python = ">=3.12"
dependencies = [
    "requests>=2.31.0",
]

# Standard PEP 621 extras definition
[project.optional-dependencies]
gui = ["PyQt6"]

# uv manages lockfile settings here
[tool.uv]

```

### Scripts Handling 

#### 1. Defining the Script (Config)

You stop using `[tool.poetry.scripts]` and switch to the standard **PEP 621** format in `pyproject.toml`.

**Poetry:**

```toml
[tool.poetry.scripts]
myapp = "my_package.main:run"

```

**uv (Standard):**

```toml
[project.scripts]
myapp = "my_package.main:run"

```

#### 2. Installing the Script

When you run `uv sync`, uv installs your current project into the virtual environment (in editable mode by default).

* **Command:** `uv sync` (or `uv sync --all-extras`)
* **Result:** The `myapp` executable is generated and placed inside `.venv/bin/` (Mac/Linux) or `.venv\Scripts\` (Windows).

#### 3. Running the Script

You have two ways to run it, equivalent to Poetry:

##### Method A: The Direct Equivalent (`uv run`)

This works exactly like `poetry run`. It ensures the environment is set up and then executes the command.

```bash
uv run myapp

```

* **⚠️ Note on Extras:** As discussed, if your script relies on "extras" dependencies, `uv run` is strict. It will strip extras if you don't explicitly ask for them.
* *If you need extras:* `uv run --all-extras myapp`



##### Method B: The "Activated" Way (Shell)

Since `uv sync` installs the script into the `.venv` `bin` folder, you can just run it directly if your environment is active.

```bash
# 1. Activate
source .venv/bin/activate

# 2. Run directly (No 'uv' prefix needed)
myapp

```

* **Benefit:** This **preserves your installed extras**. If you ran `uv sync --all-extras` previously, running `myapp` directly respects that state without you typing flags every time.

#### Summary Cheat Sheet

| Action             | Poetry                    | uv                                     |
|--------------------|---------------------------|----------------------------------------|
| **Define Script**  | `[tool.poetry.scripts]`   | `[project.scripts]`                    |
| **Install Script** | `poetry install`          | `uv sync`                              |
| **Run (Wrapper)**  | `poetry run myapp`        | `uv run myapp`                         |
| **Run (Shell)**    | `poetry shell` -> `myapp` | `source .venv/bin/activate` -> `myapp` |

---

### 🛠 Workflow: Setting up a new project

1. **Initialize:**
```bash
uv init my-project
cd my-project

```


2. **Create Venv (Optional, `uv sync` does this automatically):**
```bash
uv venv

```


3. **Add Dependencies:**
```bash
uv add fastAPI
uv add --dev ruff
# Add dependency to an 'extra' group named 'gui'
uv add PyQt6 --optional gui

```


4. **Sync (Install with Extras):**
```bash
# IMPORTANT: Include extras here, otherwise they are removed
uv sync --all-extras

```


5. **Run:**
```bash
uv run --all-extras main.py

```



---

### 📝 Important Notes

1. **`uv sync` is the new `poetry install**`:
`uv sync` is the most important command. It ensures your `.venv` exactly matches your `uv.lock`.
2. **⚠️ The Trap of `uv sync` and Extras**:
`uv sync` is **strictly declarative**. If you have extras installed and you run a plain `uv sync` (without `--all-extras` or `--extra ...`), **uv will uninstall your extras** to match the "base" definition. You must include the flag every time you sync if you want them present.
3. **Performance**:
uv is written in Rust and is significantly faster than Poetry. You will notice `uv add` and `uv sync` are nearly instant.
4. **Python Version Management**:
uv can install Python versions for you.
`uv python install 3.12` -> `uv venv --python 3.12`
5. **Scripts**:
Instead of `[tool.poetry.scripts]`, use `[project.scripts]`.