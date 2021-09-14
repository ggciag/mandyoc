# Contributing to MANDYOC

Contributions to MANDYOC are welcome.
If you have an issue, a bug, a code contribution or a documentation
contribution, **thanks for helping to improve MANDYOC!**

**How to contribute to MANDYOC?**

- Submitting bug reports and feature requests
- Writing tutorials or examples
- Fixing typos and improving to the documentation
- Writing code for everyone to use

## Reporting a Bug

Find the Issues and click New Issue.
Remember to _choose bug label_ for this type of issues.

When creating a bug report, please be as specific as possible when describing
how to reproduce an issue, and include both the intended/expected result and
what you are actually getting.
Attach as much as possible of the following information to your issue:

- a minimal parameter file that reproduces the issue,
- the log.txt file that was created during the model run,
- the error message you saw on your screen,
- any information that helps us understand why you think this is a bug, and
  how to reproduce it.

## Making MANDYOC Better

**First off, thank you for considering contributing to our project!**

If you want to make some kind of contribution to MANDYOC, please note the
following general guidelines.

### General guidelines

We _follow the
[git pull request workflow](https://www.asmeurer.com/git-workflow/) to make
changes to our code base_.
Every change made goes through a pull request, even our own, so that our
continuous integration services have a change to check that the code is up to
standards and passes all our tests.
This way, the master/main branch is always stable.

General guidelines for pull requests (PRs):

- **Open an issue first describing what you want to do**.
  If there is already an issue that matches your PR, leave a comment there
  instead to let us know what you plan to do.
- **Create a fork** of the repository on your personal account.
  **Then, clone it** in your computer:

  ```bash
  git clone your-fork-url
  ```

- **Create a separate branch** (sometimes called a feature branch) on which
  you do your modifications:

  ```bash
  git checkout -b branch-name
  ```

- Once you have created your branch, **make your changes and commit them**.
  Remember to keep your commits atomic, that is, each commit should represent a
  single unit of change.
  Also, remember to write helpful commit messages, so that someone can
  understand what the commit does.

  ```bash
  git add filename

  git commit
  ```
  After that **push up your changes** to your fork:

  ```bash
  git push
  ```

- **Make a pull request**.
  Ensure the PR description clearly describes the problem and solution.
  Include the relevant issue number.

**Remember to test your changes**. We have some global test code to check that 
the Mandyoc output result for the Crameri model is equal to the expected result.
To make this test you need to install `pytest` and run:

```bash
make test
```

## Do you have questions about MANDYOC code?

Ask any question you have by **opening up an Issue**, and labeling it as
a _question_.
