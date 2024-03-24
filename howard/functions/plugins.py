import importlib
import os


def load_plugins(plugins_dir: str) -> dict:
    """
    The function `load_plugins` loads Python plugins from a specified directory and returns a dictionary
    mapping plugin names to their respective modules.

    :param plugins_dir: The `plugins_dir` parameter in the `load_plugins` function is a string that
    represents the directory where the plugins are located. This function loads Python plugins from the
    specified directory and returns a dictionary containing the loaded plugins
    :type plugins_dir: str
    :return: A dictionary containing information about the loaded plugins is being returned. Each key in
    the dictionary represents the name of a plugin, and the corresponding value is a dictionary
    containing the attributes and functions defined in that plugin.
    """

    plugins = {}
    # plugin_dir = os.path.join(os.path.dirname(__file__), plugins_dir)

    for plugin_name in os.listdir(plugins_dir):
        plugin_path = os.path.join(plugins_dir, plugin_name)

        if os.path.isdir(plugin_path) and os.path.exists(
            os.path.join(plugin_path, "__init__.py")
        ):
            spec = importlib.util.spec_from_file_location(
                f"plugins.{plugin_name}", os.path.join(plugin_path, "__init__.py")
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            plugin_dict = {key: value for key, value in list(module.__dict__.items())}
            plugins[plugin_name] = plugin_dict

    return plugins


def list_plugins(plugins: dict, plugins_dir: str) -> dict:
    """
    The function `list_plugins` loads plugin information from a specified directory and checks which
    plugins are enabled.

    :param plugins: The `plugins` parameter is a dictionary containing information about various
    plugins. Each key in the dictionary represents the name of a plugin, and the corresponding value is
    a dictionary containing data about that plugin
    :type plugins: dict
    :param plugins_dir: The `plugins_dir` parameter is a string that represents the directory where the
    plugins are located. In the `list_plugins` function, this directory is used to locate the plugins
    and gather information about them
    :type plugins_dir: str
    :return: The function `list_plugins` returns a dictionary containing information about each plugin
    specified in the `plugins` parameter. The information includes whether the plugin is enabled (based
    on whether it is in the list of enabled plugins obtained from the specified directory), as well as
    any additional data provided for each plugin in the `plugins` dictionary.
    """

    plugin_info = {}
    # plugin_dir = os.path.join(os.path.dirname(__file__), plugins_dir)
    plugins_init = os.path.join(plugins_dir, "__init__.py")
    spec = importlib.util.spec_from_file_location(f"plugins", plugins_init)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    plugins_enabled = getattr(module, "__all__", [])
    for plugin_name, plugin_data in plugins.items():
        plugin_info[plugin_name] = {
            "enabled": plugin_name in plugins_enabled,
            **plugin_data,
        }

    return plugin_info


def plugins_to_load(list_plugins_dict: dict) -> dict:
    """
    The function `plugins_to_load` filters a dictionary of plugins based on their "enabled" and
    "__enabled__" keys.

    :param list_plugins_dict: The `list_plugins_dict` parameter is a dictionary containing information
    about various plugins. Each key in the dictionary represents the name of a plugin, and the
    corresponding value is another dictionary with plugin information. The plugin information dictionary
    may contain keys such as "enabled" and "__enabled__" to indicate whether the
    :type list_plugins_dict: dict
    :return: The function `plugins_to_load` returns a dictionary containing plugins that are enabled
    based on the input `list_plugins_dict`. The plugins are selected based on the values of the
    "enabled" and "__enabled__" keys in the nested dictionaries within the input dictionary.
    """
    plugins_to_load = {}
    for plugin_name in list_plugins_dict:
        if list_plugins_dict[plugin_name].get("enabled", False) and list_plugins_dict[
            plugin_name
        ].get("__enabled__", False):
            plugins_to_load[plugin_name] = list_plugins_dict[plugin_name]

    return plugins_to_load
