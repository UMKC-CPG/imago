---
layout: default
title: Script Catalog
order: 99
has_toc: false
---

# Script Catalog

The Script catalog is simply a list of all the available scripts that come with Imago and a simple description of their purpose. This is meant to be searchable and the descriptions are designed to include key words about functionality.

## Catalog:

{% for page in site.reference_manual %}
   {% if page.script %}
   <h3><a href="{{ page.url }}"> {{ page.script }} </a></h3>

   {{ page.description }}
   {% endif %}
{% endfor %}
