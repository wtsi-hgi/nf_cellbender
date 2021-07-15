#!/usr/bin/env bash

git add * && \
  git add .gitlab-ci.yml && \
  git add .gitignore && \
  git commit --allow-empty -m "auto_gitlab; ci_run_kill; ci_run_lsf" && \
  git push -u origin
