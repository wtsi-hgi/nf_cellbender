#!/usr/bin/env bash

git add * && \
    git add .gitignore && \
    git add .gitlab-ci.yml && \
    git commit -m "auto_gitlab" && \
    git push -u origin
