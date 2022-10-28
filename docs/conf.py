# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import re

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "mvMAPIT"
copyright = "2022, lcrawlab"
author = "lcrawlab"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "m2r2",
    "sphinx_immaterial",
    "sphinx.ext.todo",
    "sphinx_jinja",
]

# Get the list of all files and directories
path = "../man"
dir_list = os.listdir(path)
files = {}
for f in dir_list:
    fp = os.path.join(path, f)
    with open(fp, "r") as rd:
        doc = rd.read()

        doctype_search = re.search("docType{(.*)}", doc, re.IGNORECASE)
        if doctype_search:
            continue

        title_search = re.search("title{([^}]*)}", doc, re.IGNORECASE)
        if title_search:
            title = title_search.group(1)
        else:
            title = f

        name_search = re.search("name{(.*)}", doc, re.IGNORECASE)
        if name_search:
            name = name_search.group(1)
        else:
            name = f

        parameter_search = re.findall(r"item{([^}]*)}{([^}]*)}", doc)
        if parameter_search:
            parameters = {}
            parameters = dict(parameter_search)
            for key, value in parameters.items():
                parameters[key] = value.replace("\n", " ")

        usage_search = re.search(r"usage{([^}]*)}", doc, re.MULTILINE)
        if usage_search:
            usage = usage_search.group(1)
        else:
            continue

        value_search = re.search(r"value{([^}]*)}", doc, re.MULTILINE)
        if value_search:
            value = value_search.group(1)
        else:
            continue

        description_search = re.search(r"description{([^}]*)}", doc, re.MULTILINE)
        if description_search:
            description = description_search.group(1)
        else:
            continue

        files[f] = {
            "name": name,
            "usage": usage,
            "name_len": len(name),
            "parameters": parameters,
            "title": title,
            "value": value,
            "description": description,
        }

jinja_contexts = {
    "file_list": {"rd": files},
}


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_immaterial"  # "alabaster"
html_static_path = ["_static"]

# The name for this set of Sphinx documents.  If None, it defaults to
html_static_path = ["_static"]
# html_css_files = ["extra_css.css"]
html_last_updated_fmt = ""
html_title = "mvMAPIT"
html_favicon = "_static/images/favicon.ico"
html_logo = "_static/images/favicon-w.png"

html_theme_options = {
    "icon": {
        "repo": "fontawesome/brands/github",
    },
    "site_url": "https://lcrawlab.github.io/mvMAPIT/",
    "repo_url": "https://github.com/lcrawlab/mvMAPIT/",
    "repo_name": "mvMAPIT",
    "repo_type": "github",
    "edit_uri": "blob/main/docs",
    "globaltoc_collapse": True,
    "features": [
        "navigation.expand",
        # "navigation.tabs",
        # "toc.integrate",
        "navigation.sections",
        # "navigation.instant",
        # "header.autohide",
        "navigation.top",
        # "navigation.tracking",
        # "search.highlight",
        "search.share",
        "toc.follow",
        "toc.sticky",
    ],
    "palette": [
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "indigo",
            "accent": "green",
            "toggle": {
                "icon": "material/lightbulb-outline",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "grey",
            "accent": "green",
            "toggle": {
                "icon": "material/lightbulb",
                "name": "Switch to light mode",
            },
        },
    ],
    # BEGIN: version_dropdown
    "version_dropdown": True,
    "version_info": [
        {
            "version": "https://lcrawlab.github.io/mvMAPIT",
            "title": "Github Pages",
            "aliases": [],
        },
    ],
    # END: version_dropdown
    "toc_title_is_page_title": True,
    # BEGIN: social icons
    "social": [
        {
            "icon": "fontawesome/brands/github",
            "link": "https://github.com/lcrawlab/mvMAPIT",
        },
        # {
        #    "icon": "fontawesome/brands/python",
        #    "link": "https://pypi.org/project/sphinx-immaterial/",
        # },
    ],
    # END: social icons
}

# Output file base name for HTML help builder.
htmlhelp_basename = "mvMAPIT-doc"
