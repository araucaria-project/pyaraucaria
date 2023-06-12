from lark import Lark, Tree
from lark.exceptions import UnexpectedCharacters
import logging
import sys


logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__.rsplit('.')[-1])


class ObsPlanParser:

    @property
    def _line_grammar(self):
        lin_gr = r"""
        ?start          : sequences
        !sequences      : sequence*
        !sequence       : end_line* begin_sequence args kwargs comment? separator all_commands end_sequence comment? (end_line* comment?)*
        !all_commands   : (command separator)*
        !command        : (command_name args kwargs | sequence) comment?
        !command_name   : word
        kwargs          : kwarg*
        args            : val*
        kwarg           : kw "=" (val ",")+ val | kw "=" val
        !begin_sequence : "BEGINSEQUENCE"
        !end_sequence   : "ENDSEQUENCE"
        separator       : end_line+
        word            : /[^{BEGINSEQUENCE}{ }][A-Z]+/
        comment         : /#.*/
        end_line        : /\n/
        !kw             : string_simple
        !val            : string_simple | string_quoted
        ?string_quoted  : /[\"].*[\"]|[\'].*[\']/
        ?string_simple  : /[^\s=\'\"]+/
        %ignore /[ \f]+/
        %ignore "#" /[^\n]/*
        """
        return lin_gr

    @staticmethod
    def _build_kwargs(tree: Tree):

        kwargs_dict = {}
        for child in tree.children:
            if child.data == 'kwarg':
                kw = ''
                val = ''
                for child_2 in child.children:
                    if child_2.data == 'kw':
                        kw = str(child_2.children[0])
                    if child_2.data == 'val':
                        val = str(child_2.children[0])
                kwargs_dict[kw] = val
        return kwargs_dict

    @staticmethod
    def _build_args(tree: Tree):

        all_args_list = []
        for child in tree.children:
            all_args_list.append(str(child.children[0]))
        return all_args_list

    @staticmethod
    def _build_command_name(tree: Tree) -> str:

        for child in tree.children:
            if child.data == 'word':
                word = str(child.children[0])
                return word

    def _build_command(self, tree: Tree):

        command_dict = {}
        for child in tree.children:
            if child.data == 'command_name':
                com_name = self._build_command_name(child)
                if com_name == 'sKYFLAT':
                    com_name = 'SKYFLAT'
                if com_name == 'sTOP':
                    com_name = 'STOP'
                command_dict['command_name'] = com_name
            if child.data == 'args':
                if child.children:
                    if child.children[0].data == 'val':
                        command_dict['args'] = self._build_args(child)
            if child.data == 'kwargs':
                if child.children:
                    if child.children[0].children[1].data == 'val':
                        command_dict['kwargs'] = self._build_kwargs(child)
        return command_dict

    def _build_all_commands(self, tree: Tree):

        all_commands_list = []
        for child in tree.children:
            if child.data == 'command':
                if child.children[0].data == 'command_name':
                    all_commands_list.append(self._build_command(child))
                if child.children[0].data == 'sequence':
                    all_commands_list.append(self._build_sequence(child.children[0]))
        return all_commands_list

    def _build_sequence(self, tree: Tree):

        sequence_dict = {}
        for child in tree.children:
            if child.data == 'begin_sequence':
                sequence_dict['begin_sequence'] = 'begin'
            if child.data == 'args':
                if child.children:
                    if child.children[0].data == 'val':
                        sequence_dict['args'] = self._build_args(child)
            if child.data == 'kwargs':
                if child.children:
                    if child.children[0].children[1].data == 'val':
                        sequence_dict['kwargs'] = self._build_kwargs(child)
            if child.data == 'all_commands':
                sequence_dict['all_commands'] = self._build_all_commands(child)
        return sequence_dict

    def _build_sequences(self, tree: Tree):

        sequence_list = []
        try:
            for child in tree.children:
                if child.data == 'sequence':
                    sequence_list.append(self._build_sequence(child))
            return sequence_list
        except AttributeError:
            log.error(f'Text cannot be parsed, please check string')

    def _parse_text(self, text: str) -> Tree:

        line_parser = Lark(self._line_grammar)
        parse = line_parser.parse
        try:
            parsed = parse(self._prepare_text(text))
            return parsed
        except AttributeError:
            log.error(f'Text cannot be parsed, please check string')
        except UnexpectedCharacters:
            log.error(f'No terminal matches in the current parser context')

    @staticmethod
    def _prepare_text(text: str) -> str:
        """
        Prepare text to parse. If BEGINSEQUENCE and ENDSEQUENCE is not added it will add.
        Also solve SKYFLAT lark error and replace for sKYFLAT, later replace back.
        :param text: text to parse
        :return: prepared text to parse
        """
        return_text = text
        log.debug(f'{return_text.find("SKYFLAT")}')
        if text.find('BEGINSEQUENCE') < 0:
            return_text = f"BEGINSEQUENCE \n{text} \nENDSEQUENCE"
        if return_text.find('SKYFLAT') >= 0:
            return_text = return_text.replace('SKYFLAT', 'sKYFLAT')
        if return_text.find('STOP') >= 0:
            return_text = return_text.replace('STOP', 'sTOP')
        return return_text

    @staticmethod
    def _read_file(file_name: str) -> str:
        """
        Read txt file
        :param file_name: file name
        :return: string
        """
        try:
            file = str(open(file_name, "r").read())
            return file
        except FileNotFoundError:
            log.error(f'File: {file_name} not found.')

    @staticmethod
    def _write_to_file(file_name: str, sequences) -> None:

        file = open(file_name, "w")
        log.debug(f'{sequences}')
        file.write(str(sequences))
        file.close()

    def convert_from_string(self, input_string: str):
        """
        Method convert string to observation plan in json format
        :param input_string: input string
        :return: return parsed observation plan in json format
        """
        par_txt = self._parse_text(input_string)
        sequences = self._build_sequences(par_txt)
        return sequences

    def convert_from_file(self, input_file_name: str, output_file_name: str) -> None:
        """
        Method parse observation plan from txt file to txt file.
        :param input_file_name: input file name
        :param output_file_name: output file name
        :return:
        """
        text = self._read_file(input_file_name)
        par_txt = self._parse_text(text)
        sequences = self._build_sequences(par_txt)
        self._write_to_file(output_file_name, sequences)


if __name__ == '__main__':
    opp = ObsPlanParser()
    try:
        txt = sys.argv[1]
        log.debug(txt)
        result = opp.convert_from_string(txt)
        if result:
            log.info(f'{result}')
    except IndexError:
        log.error(f'Please add string to conversion')
