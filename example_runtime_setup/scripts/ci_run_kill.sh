#!/usr/bin/env bash

git add * && \
  git add .gitlab-ci.yml && \
  git add .gitignore && \
  git commit --allow-empty -m "auto_gitlab; ci_run_kill" && \
  git push -u origin
