import importlib
import os


def plugins_infos(plugins_dir: str, subfolder_plugins: str = "plugins") -> dict:
    """
    The `plugins_infos` function loads Python plugins from a specified directory and returns a
    dictionary mapping plugin names to their respective modules.

    :param plugins_dir: The `plugins_dir` parameter in the `plugins_infos` function is a string that
    represents the directory where the plugins are located. This function loads Python plugins from the
    specified directory and returns a dictionary containing the loaded plugins
    :type plugins_dir: str
    :param subfolder_plugins: The `subfolder_plugins` parameter in the `plugins_infos` function is a
    string that represents the subfolder within the `plugins_dir` where the plugins are located. By
    default, the value of `subfolder_plugins` is set to "plugins". This parameter is used to specify the
    subfolder, defaults to plugins
    :type subfolder_plugins: str (optional)
    :return: A dictionary containing information about the loaded plugins is being returned. Each key in
    the dictionary represents the name of a plugin, and the corresponding value is a dictionary
    containing the attributes and functions defined in that plugin.
    """

    plugins = {}

    # For each plugin folder
    for plugin_name in os.listdir(plugins_dir):

        # plugin path
        plugin_path = os.path.join(plugins_dir, plugin_name)

        # If plugin is dir and __init__.py exists
        if os.path.isdir(plugin_path) and os.path.exists(
            os.path.join(plugin_path, "__init__.py")
        ):
            # plugin spec
            spec = importlib.util.spec_from_file_location(
                f"{subfolder_plugins}.{plugin_name}",
                os.path.join(plugin_path, "__init__.py"),
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)

            # Create plugin info dict
            plugin_dict = {key: value for key, value in list(module.__dict__.items())}
            plugins[plugin_name] = plugin_dict

    return plugins


def plugins_list(
    plugins: dict, plugins_dir: str, subfolder_plugins: str = "plugins"
) -> dict:
    """
    The `plugins_list` function loads plugin information from a specified directory and determines which
    plugins are enabled based on a dictionary of plugin data.

    :param plugins: The `plugins` parameter is a dictionary containing information about various
    plugins. Each key in the dictionary represents the name of a plugin, and the corresponding value is
    a dictionary containing data about that plugin
    :type plugins: dict
    :param plugins_dir: The `plugins_dir` parameter is a string that represents the directory where the
    plugins are located. This directory is used by the `list_plugins` function to locate the plugins and
    gather information about them
    :type plugins_dir: str
    :param subfolder_plugins: The `subfolder_plugins` parameter in the `plugins_list` function is a
    string that represents the subfolder within the `plugins_dir` where the plugins are located. By
    default, the value of `subfolder_plugins` is set to "plugins". This parameter is used to specify the
    subfolder, defaults to plugins
    :type subfolder_plugins: str (optional)
    :return: The function `plugins_list` returns a dictionary `plugin_info` containing information about
    each plugin specified in the `plugins` parameter. The information includes whether the plugin is
    enabled (based on whether it is in the list of enabled plugins obtained from the specified
    directory), as well as any additional data provided for each plugin in the `plugins` dictionary.
    """

    plugin_info = {}

    # Plugins init
    plugins_init = os.path.join(plugins_dir, "__init__.py")

    if os.path.exists(plugins_init):

        # Module plugins spec
        spec = importlib.util.spec_from_file_location(subfolder_plugins, plugins_init)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # List of enabled pugins in module
        plugins_enabled = getattr(module, "__all__", [])

        # Add enable for each plugins
        for plugin_name, plugin_data in plugins.items():
            plugin_info[plugin_name] = {
                "enabled": plugin_name in plugins_enabled,
                **plugin_data,
            }

    return plugin_info


def plugins_to_load(plugins_list_dict: dict) -> dict:
    """
    The `plugins_to_load` function filters a dictionary of plugins based on their "enabled" and
    "__enabled__" keys.

    :param plugins_list_dict: The `plugins_list_dict` parameter is a dictionary containing information
    about various plugins. Each key in the dictionary represents the name of a plugin, and the
    corresponding value is another dictionary with plugin information. The plugin information dictionary
    may contain keys such as "enabled" and "__enabled__" to indicate whether the
    :type plugins_list_dict: dict
    :return: The function `plugins_to_load` returns a dictionary containing plugins that are enabled
    based on the input `plugins_list_dict`. The plugins are selected based on the values of the
    "enabled" and "__enabled__" keys in the nested dictionaries within the input dictionary.
    """

    plugins_to_load = {}
    for plugin_name in plugins_list_dict:
        if plugins_list_dict[plugin_name].get("enabled", False) and plugins_list_dict[
            plugin_name
        ].get("__enabled__", False):
            plugins_to_load[plugin_name] = plugins_list_dict[plugin_name]

    return plugins_to_load
