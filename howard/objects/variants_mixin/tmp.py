from howard.functions.commons import (
    get_tmp,
    remove_if_exists
)

class variants_tmp:

    def set_tmp_files(self, tmp_files: list = None) -> None:
        """
        Set temporary files

        :param tmp_files: The `tmp_files` parameter in the `set_tmp_files` method is a list of temporary
        files that you want to set as the `tmp_files` attribute of the class. If no files are provided,
        it defaults to an empty list
        :type tmp_files: list

        :return: The `set_tmp_files` method is returning `None`.
        """

        # Init
        if tmp_files is None:
            tmp_files = []

        # Create tmp_files attribute
        self.tmp_files = tmp_files

    def get_tmp_files(self) -> list:
        """
        Get temporary files
        """

        # Return tmp_files attribute
        return self.tmp_files

    def add_tmp_files(self, tmp_files: list = None) -> None:
        """
        Extend a temporary file list to the list of temporary files

        :param tmp_files: The `tmp_files` parameter in the `add_tmp_files` method is a list of temporary
        files that you want to add to the `tmp_files` attribute of the class. If no files are provided,
        it defaults to an empty list
        :type tmp_files: list

        :return: The `add_tmp_files` method is returning `None`.
        """

        # Init
        if tmp_files is None:
            tmp_files = []

        # Append list to tmp_files attribute
        self.tmp_files.extend(tmp_files)

    def remove_tmp_files(self, tmp_files: list = None) -> None:
        """
        Remove files from temporary files attribute

        :param tmp_files: The `tmp_files` parameter in the `remove_tmp_files` method is a list of
        temporary files that you want to remove from the `tmp_files` attribute of the class. If no
        files are provided, it defaults to an empty list
        :type tmp_files: list

        :return: The `remove_tmp_files` method is returning `None`.
        """

        # Init
        if tmp_files is None:
            tmp_files = []

        # Remove tmp files from tmp_files attribute
        self.tmp_files = [f for f in self.get_tmp_files() if f not in tmp_files]

    def clean_tmp_files(self) -> list:
        """
        Remove temporary files
        """

        # Remove tmp_files files
        tmp_files = self.get_tmp_files()
        remove_if_exists(tmp_files)
        self.tmp_files = []

        return tmp_files

    def get_tmp_dir(self) -> str:
        """
        The function `get_tmp_dir` returns the temporary directory path based on configuration
        parameters or a default path.
        :return: The `get_tmp_dir` method is returning the temporary directory path based on the
        configuration, parameters, and a default value of "/tmp".
        """

        return get_tmp(
            config=self.get_config(), param=self.get_param(), default_tmp="/tmp"
        )
