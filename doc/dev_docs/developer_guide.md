# Developer guide

## Introduction
The goal of this document is to describe requirements and best practices for contributing to the RAMSES codebase.
It is intended for developers and contributors who are familiar with the codebase and want to make changes or add new features.
The scope of this guide is technical, and it does not cover the scientific aspects of the code.

## Branching model summary

RAMSES uses two main branches: `dev` and `stable`. `dev` is the development branch where new features are added before being put into production on the `stable` branch.
The `dev` branch is meant for developers, while the `stable` version is meant for users.

**All** new features, bug fixes, and improvements should stem from and be merged into the `dev` version.
When starting to work on a new project, we encourage the following procedure, which will be detailed further in the next sections:
1. Create a private fork of the `ramses-organisation/ramses` repository  (section 3.1).
2. Create a new branch of `dev` with a descriptive name  (section 3.1).
3. Make your changes and push them to your fork (section 3.2).
4. Once ready, submit a pull request to the `ramses-organisation/ramses` repository that points onto `dev` (section 3.3).


Periodically, the developments on the `dev` branch that have proven to be robust, will be merged into the `stable` branch (see section 3.4). What is merged and what not, is decided through collegial discussions.
At any point, any commit in `stable` should be a subset of the commits in `dev`.

## Development workflow

### Creating your working branch

First, developers create their own private fork of the main repository. More information on forks and how to create them on github can be found [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo).

We provide here some example commands to implement this workflow:
```
# Clone the public version
git clone https://github.com/ramses-organisation/ramses.git
# Move into the directory
cd ramses
# Assuming you have already forked the repository, add it as a remote
git remote add fork git@github.com:your_username/ramses.git
# Switch to the dev branch
git checkout dev
# Create a new branch
git checkout -b my-new-feature
```

### Committing changes

```
# Switch to your working branch
git checkout my-new-feature
# Make your changes
git add [path to updated files]
# Validate your changes
git commit -m "Add my new feature"
# Push the changes to your fork (and optionally set it as the upstream)
git push --set-upstream fork my-new-feature

# Keep in sync with the public dev branch
# This should be done on a regular basis
git pull --rebase --autostash origin dev

# Make more changes
git add [path to updated files]
git commit -m "Add more features"
git push  # since we have set the upstream, this will push onto your fork
```

We do not enforce a strict commit message format, but we encourage the following guidelines:
- **Short summary**: A one-line summary of the changes.
- **Detailed description**: A more detailed description of the changes.
- **Commit frequency**: We encourage frequent commits, but not too frequent. A good rule of thumb is to commit whenever a logical unit of work is completed and the code compiles and runs without errors.


### Submitting a pull request

Once you are ready to submit your changes, you can create a pull request (PR) on the GitHub interface.
Any new feature, bug fix, or improvement should be submitted as a PR to the `dev` branch.
A PR should contain one logical unit of work. If you have multiple features or bug fixes, please submit them as separate PRs.
More information on how to submit a PR on github, see [the github doc](
https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

In order to be merged, a PR
- should be reviewed by at least one developer,
- should pass the continuous integration tests,
- have some information about the changes in the description,
- include documentation (for new features),
- include tests (for new features)
- should be in sync with the current state of the `dev` branch (either va `git pull --rebase --autostash origin dev` or using the appropriate button in the Github interface).

It is the responsibility of the reviewer to ensure these conditions are met before merging the PR.

We encourage the use of the GitHub interface to tag the PR with the appropriate labels (e.g., `enhancement`, `bug`, `on-merge: backport-to-stable`). The tag `on-merge: backport-to-stable` will automatically trigger backporting of this PR onto the `stable` branch once merged. We discourage manual backporting of bugfixes.

### Updates to `stable`

No pull-request nor commit should be pushed onto the `stable` branch, with the exception of:
- (automated) backported commits,
- regular updates of the `stable` version.
In order to enforce this, the branch is protected (no one should commit directly onto it).

**Backported commits** include *small* bugfixes and modification that should be included in the stable version (e.g., update of the contributor list). It excludes any breaking change, new feature, or experimental code.

**Updates of the stable version** happen on regular occasions (e.g., every Ramses User Meeting). The `stable` branch will be updated with the latest changes from `dev`. Collegial discussions will be held to decide which features are ready to be included in the `stable` version and which ones should be left out.


## GitHub Management

All developments related to the RAMSES community are hosted under the [Ramses organisation](https://github.com/ramses-organisation).
The owner of the organisation is Romain Teyssier, but the management of the organisation and the repositories is done by the RAMSES community.

**In order to encourage contributions and to promote a more inclusive and community-driven development, contributors to the code will be listed online** (to reward their effort) and be granted permission based on their level of engagement.
We distinguish three user roles:
- **Administrators**: Have full access to all repositories and can manage the organisation (`admin` role on GitHub). Their number should be kept to a minimum and they should be trusted members of the community. They are in charge of administrating the GitHub organisation. This includes notably creating/deleting repositories, managing user permissions, configuring default branches and protections on those.
- **Maintainers**: Have write access to the repositories and can merge pull requests (`maintain` role on GitHub). They are responsible for the day-to-day management of the repositories. This includes notably merging pull requests and closing solved issues. They should be active members of the community with a good understanding of the codebase.
- **Contributors**: Can manage issues and pull requests (`triage` role on GitHub). They are in charge of reviewing pull requests, suggest changes and approve pull requests, as well as answering issues. Anyone who has made a significant contribution to RAMSES, including code contribution or community management (e.g. creating GitHub issues or responding to other users). This role rewards active members of the community, and should be used to encourage new contributors, especially early-career researchers.

Promotion to a higher role is done by one of the administrators. The decision should be based on the activity and the quality of the contributions of the user. It is decided by a vote of the community (when promoting to contributor), by maintainers+administrators (when promoting to maintainer) or by administrators (when promoting to administrator).
Voting can be triggered by any member of the community by sending an email to the RAMSES user list, and the decision should be made within two weeks and announced on the same list. We propose a template for the email below.
To promote an individual to contributors, approval of at least two individuals is sufficient.
To promote an individual to maintainer or administrator, the approval of the voting majority is required.

Template for the email to propose to promote a user to contributor status:
```
Subject: [VOTE] Promotion of [username] to contributor

Dear RAMSES community,

I would like to propose the promotion of [username] to the role of contributor.
[username] has been an active member of the community for [time], and has made
significant contributions (contribution 1, contribution 2, ...). I believe
that [username] would be a valuable addition to the team.

Best regards,
[your name]
```

To promote a user to maintainer or administrator, the email should be sent to the RAMSES user list (ramses_users@googlegroups.com) with the subject `[VOTE] Promotion of [username] to maintainer` or `[VOTE] Promotion of [username] to administrator`, respectively. Replies should be collected for a week by one of the administrators, either through direct reply to the email or using some online platform (for example https://framadate.org/abc/en/), and the decision should be announced on the same list.


Any member of the community can demand to demote a user to a lower role. Valid reasons to demote from administrator or maintainer include:
- Absence of activity for a long period of time (e.g., no response to emails),
- Repetitive mistakes in the management of the repositories (e.g., deleting important branches, merge without review, ignoring this guide).
Demotion from the contributor list should only be done only in extreme cases (e.g., harassment, spamming) and will be discussed by the administrators on a case-by-case basis.

Finally, the list of current administrators, maintainers, and contributors should be kept up-to-date in the [`CONTRIBUTORS.md`](./contributors.md) file of the ramses-organisation repository. This list should be updated whenever a new user is promoted or demoted and include:
- Name of the user,
- Role (administrator, maintainer, contributor),
- Date of promotion,
- List of contributions (at the time of the promotion).
