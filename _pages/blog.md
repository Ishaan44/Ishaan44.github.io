---
layout: clean
title: Blog
permalink: /blog/
---

## Welcome to my blog page!

Here you will find some of my thoughts on topics in mathematics and beyond that interest me at any given time.

# Blog Posts

{% if site.posts.size > 0 %}

{% assign posts_by_month = site.posts | group_by_exp: "post", "post.date | date: '%B, %Y'" %}

{% for month in posts_by_month %}

## {{ month.name }}

{% for post in month.items %}

### [{{ post.title }}]({{ post.url | relative_url }})

{% if post.description %}
{{ post.description }}
{% endif %}

[Read more →]({{ post.url | relative_url }})

{% endfor %}

{% endfor %}

{% else %}

No blog posts yet.

{% endif %}
