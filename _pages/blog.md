---
layout: page
title: Blog
permalink: /blog/
nav: true
nav_order: 4
---

Occasional posts on probability, mathematics, and things I am trying to understand.

{% if site.posts.size > 0 %}

{% assign posts = site.posts | sort: "date" | reverse %}
{% for post in posts %}

## [{{ post.title }}]({{ post.url | relative_url }})

<span class="post-meta">{{ post.date | date: "%B %-d, %Y" }}</span>

{% if post.description %}
{{ post.description }}
{% endif %}

{% endfor %}

{% else %}

No posts yet.

{% endif %}
