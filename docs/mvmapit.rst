========
mvMAPIT
========

.. jinja:: file_list

    {% for k, v in rd.items() %}

    {{v.name}}
    {{"=" * v.name_len}}

    {{v.title | indent(4)}}

    {{v.description}}

    .. code-block:: R
        {{v.usage | indent(4)}}

    Parameters
    ----------
    {% for p, q in v.parameters.items() %}
    * **{{p}}** {{q}}
    {% endfor %}

    Return
    ------
    {{v.value}}

    ------
    {% endfor %}

Code available on github.
