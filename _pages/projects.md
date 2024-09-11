---
layout: page
title: <strong>Math</strong>
permalink: /projects/
description: Check out the details of the math I have done so far including projects, notes and other writings.
nav: true
nav_order: 3
display_categories: [Projects, Notes and other writings]
horizontal: false
---
<div class="projects">

  <!-- Projects Section -->
  <h2><strong>Projects</strong></h2>
  <ol>
    {% assign sorted_projects = site.projects | where: "category", "Projects" | sort: "importance" %}
    {% for project in sorted_projects %}
      <li>
        <div class="project-item">
          <h3>{{ project.title }}</h3>
          <p>{{ project.description | truncatewords: 20 }}</p>
          <a href="{{ project.url }}" class="read-more">Read more</a>
        </div>
      </li>
    {% endfor %}
  </ol>

  <hr>

  <!-- Research Section -->
  <h2><strong>Research</strong></h2>
  <ol>
    {% assign sorted_research = site.projects | where: "category", "Research" | sort: "importance" %}
    {% for research in sorted_research %}
      <li>
        <div class="research-item">
          <h3>{{ research.title }}</h3>
          <p>{{ research.description | truncatewords: 20 }}</p>
          <a href="{{ research.url }}" class="read-more">Read more</a>
        </div>
      </li>
    {% endfor %}
  </ol>

  <hr>

  <!-- Notes and Other Writings Section -->
  <h2><strong>Notes and Other Writings</strong></h2>
  <ol>
    {% assign sorted_notes = site.projects | where: "category", "Notes and Other Writings" | sort: "importance" %}
    {% for note in sorted_notes %}
      <li>
        <div class="note-item">
          <h3>{{ note.title }}</h3>
          <p>{{ note.description | truncatewords: 20 }}</p>
          <a href="{{ note.url }}" class="read-more">Read more</a>
        </div>
      </li>
    {% endfor %}
  </ol>

</div>
