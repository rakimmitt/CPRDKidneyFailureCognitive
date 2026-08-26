# CPRDKidneyFailureCognitive

## Project purpose

This is an academic research repository containing R code for epidemiological
analyses of cognitive outcomes in people with advanced chronic kidney disease
using CPRD Aurum and related datasets.

The repository contains analysis code, codelists, and supporting files.
Large or sensitive source datasets should not be stored in this repository.

## Working style

- The primary programming language is R.
- Prefer tidyverse-style R where this is consistent with the existing code.
- Use `<-` rather than `=` for assignment in R code.
- Prefer clear, readable code over unnecessarily compact code.
- Preserve the existing structure and naming conventions unless there is a
  good reason to change them.
- Where user-specific final database tables are required, use the `rk_`
  prefix where consistent with the surrounding code.
- Prefer relative project paths or `file.path()` where practical rather than
  introducing machine-specific absolute Windows paths.

## Before changing code

- Inspect the relevant surrounding scripts and functions before making changes.
- Do not assume what a CPRD/Aurum helper function does purely from its name;
  inspect its use in the repository where possible.
- For a non-trivial change, briefly explain the proposed approach before editing.
- Preserve existing behaviour unless the requested task explicitly changes it.
- Do not silently refactor unrelated code while making a focused change.

## R and analysis principles

- Make cohort definitions, dates, exclusions, censoring rules, and outcome
  definitions explicit.
- Be particularly careful about distinguishing prevalent from incident diagnoses.
- Do not change epidemiological definitions or statistical assumptions merely
  to make code easier to write.
- Flag situations where a proposed coding change could alter the study population,
  exposure definition, outcome definition, follow-up, or interpretation.
- When generating codelists, favour reproducible search and review steps and
  preserve enough information for manual clinical review.
- Do not remove apparently redundant clinical search terms without checking
  whether doing so could change codelist sensitivity.

## Database safety

- Treat the CPRD/Aurum database as sensitive research data.
- Never expose credentials, passwords, tokens, connection strings, or patient-level
  data in generated documentation, commits, prompts, or logs.
- Never add credentials or local configuration secrets to Git.
- Do not run destructive SQL or database operations such as DROP, DELETE, UPDATE,
  or overwriting an existing table unless explicitly requested.
- Before any operation that could modify persistent database content, state what
  will change.
- Prefer read-only inspection when investigating an unfamiliar table or workflow.

## Git and GitHub

- Do not commit, push, merge, rebase, force-push, or delete branches unless
  explicitly asked.
- It is fine to inspect `git status`, `git diff`, and Git history when helpful.
- Before a requested commit, summarise the files changed.
- Never commit secrets, credentials, patient-level data, or large raw datasets.

## Explaining code

The user is using this repository for both research and learning.

When asked what code does:
- explain the logic in plain English;
- identify the important objects created or changed;
- explain unfamiliar R syntax where relevant;
- distinguish necessary code from stylistic choices;
- point out potential errors or unintended consequences;
- do not rewrite working code unnecessarily.

When proposing replacement code, explain the important differences from the
existing version.

## Uncertainty

Do not invent the structure of datasets, tables, package functions, or variables.
If the answer can be established by inspecting the repository, inspect it.
If it cannot be established from available code, say what remains uncertain.