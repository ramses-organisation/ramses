.. raw:: html

   <!-- <div style="text-align:justify">  <!-- am i a maniac ? -->

.. warning::

   **NB:** You will need a **GitHub account** to follow this tutorial or
   contribute to RAMSES. If you don’t already have one, please sign up
   at https://github.com/. We also strongly recommend that you **setup
   ssh keys** following the GitHub documentation (go to your account
   ``settings``, into the ``SSH and GPG keys`` section, and follow
   instructions from there).

How to use Git
============

**RAMSES uses Git to manage code** development and version control. Git
is a distributed version control system designed to track changes in
source code during software development. It allows multiple developers
to collaborate efficiently, manage different versions of a project, and
revert to previous states if needed. With Git, one can work offline,
create new branches for new developments, and merge them back to the
main code seamlessly. While Git is a very powerful tool for community
development, using it in practice does require a well-defined community
strategy, and some basic understanding of Git’s inner workings. We have
summarized `the RAMSES Git strategy in a reference
document <https://ramses-organisation.readthedocs.io/en/dev/dev_docs/developer_guide.html>`__,
that we expand here with practical exercises.

There are many Git documentations online. We recommend this `simple
introduction <https://github.com/git-guides#learning-git-basics>`__ or
`that one <https://rogerdudler.github.io/git-guide/>`__, as well as our
`recorded RAMSES SNO webinar on Git Fundamentals by Noé
Brucy <https://ramses.cnrs.fr/a-webinar-on-git-fundamentals/>`__.

RAMSES Git model summary
------------------------

The public version of RAMSES is maintained by the
``ramses-organisation`` account on GitHub, and it is accessible here :
https://github.com/ramses-organisation/ramses.

The Git model of RAMSES relies on **one single branch**, called **dev**,
which contains the latest public version of the code. Everyone can
submit code updates to this branch through pull requests. These will
then go through a validation process before they are merged into the
main code or rejected. For obvious reasons, only a few selected persons
can edit directly the code on the ramses-organisation repository, create
branches there, or validate pull requests onto the public version. We
thus encourage the following model:

1. Create your own fork of the ramses-organisation/ramses repository.
2. In your fork, create a branch out of the ``dev`` branch to start a
   new development. Make changes and commit regularly.
3. While developing, make sure to stay in sync it with the upstream
   version.
4. Once ready, submit a pull request to the upstream repository.

All these steps are explained in detail below.

1. Creating a new fork
----------------------

The first thing to do is to **create your own fork of the ramses
repository**. A fork is basically a copy of a repository. The original
repository is called the **upstream** repository. In your fork, you can
experiment and make changes without affecting the upstream version. The
fork system allows you to synchronise your copy with the upstream repo,
and to push your developments upstream through pull requests once they
are done. An insightful documentation on forks is accessible `at this
link <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`__.

One can fork any repository, and **in this tutorial we will fork the
repository ``Mentuhotep-II/ramses``**, which is itself a fork of the
repository ``ramses-organisation/ramses`` that can be viewed as *the
source*.

.. container:: info

   **Exercise:** Create a fork of the RAMSES repository. Follow the
   `GitHub
   documentation <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo?tool=webui#forking-a-repository>`__
   but navigate to the forked ramses repository
   (`Mentuhotep-II/ramses <https://github.com/Mentuhotep-II/ramses>`__)
   instead of their “octocat/Spoon-Knife” example.

Once you have created your fork ``YOUR-USERNAME/ramses``, you can get
the code on your computer the usual way:

::

   # using SSH :
   git clone git@github.com:YOUR-USERNAME/ramses.git

2. Making changes in a new branch
---------------------------------

Now that you have your own copy of the ramses repository, you can start
making changes. In order to do this, a good practice is to create a new
branch from the ``dev`` branch, in your fork. This new branch should
have a name that refers to the development you are undertaking. In the
terminal, this is done as follows:

::

   cd ramses                      # go into your fork's directory
   git checkout dev               # go into the branch dev (of your fork)
   git checkout -b my-new-feature # create a new branch called my-new-feature.

At this point, the new branch ``my-new-feature`` exists only locally on
your computer. In order to make it visible/accessible online (in your
fork), you need to *push* with the following instruction:

::

   git push --set-upstream origin my-new-feature

Once this is done, you can check on GitHub page of your ramses fork
(https://github.com/YOUR-USERNAME/ramses) that the branch
``my-new-feature`` appears in the branches.

As you develop, you should regularly *commit* your changes. This helps
keeping track of your developments. With Git, a *commit* affects the
history on your local computer only, and you also need to *push*
developments to the original repository (i.e. your fork) online. We
provide a list of basic Git instructions below that you will need.

.. warning::

   **Basic Git instructions**

   **``git branch``** tells you on which branch you are and lists the
   local branches (use ``git branch -r`` to see the list of remote
   branches).

   **``git switch branch-name``** goes to branch ``branch-name`` (you
   can also use **``git checkout branch-name``** for that).

   **``git checkout -b new-branch``** creates a new branch (from the one
   you were in).

   **``git status``** lists the files that have been modified in the
   repository. This list also tells you which files are *staged for
   commit* (i.e. will be included in the next commit), and which are
   not.

   **``git add path/filename``** includes the file ``path/filename`` in
   the next commit.

   **``git commit -m "message"``** commits changes (staged for commit
   with ``git add``) to your local history. The message passed in quotes
   should provide a quick description of the code modifications.

   **``git push``** pushes your local history to the original online
   repository.

.. container:: info

   **Exercise:** Create a new branch. On that branch, add a comment in
   the file ``amr/constants.f90``, commit and push your changes. Use the
   GitHub web interface to check that your changes have been saved.
   :::spoiler Solution After editing the file and saving your changes
   use the following commands.

   ::

      git status    # check which files have been modified
      git add amr/constants.f90    # add amr/constants.f90 to the next commit.
      git commit -m "Added a test comment to amr/constants.f90."  # make the commit.
      git push      # push the commit online.

3. Staying in sync with the public version
------------------------------------------

The public version of the code is updated regularly, and it may happen
that it changes while you are making changes on your fork. It is
important to sync the ``dev`` branch of your fork regularly with the
original (*upstream*) version. This will make things much easier later
if you want to share your developments with the community, and it is
also necessary to benefit from bug fixes.

3.1. Synching your fork
~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to sync your fork with the upstream version is to click
on the “Sync fork” button on the GitHub page of your fork.

You can also do this in the terminal as follows. First, setup the
upstream version as the public version of RAMSES.

::

   cd ramses       # go to your local repository
   git remote add upstream git@github.com:ramses-organisation/ramses.git
   git remote -v

This should print the following lines on the terminal:

::

   origin git@github.com:YOUR-USERNAME/ramses.git (fetch)
   origin git@github.com:YOUR-USERNAME/ramses.git (push)
   upstream git@github.com:ramses-organisation/ramses.git (fetch)
   upstream git@github.com:ramses-organisation/ramses.git (push)

Once this is done, you can sync with the following command.

::

   git fetch upstream    # fetch the latest changes
   git checkout dev      # go to the dev branch
   git pull --rebase --autostash upstream dev   # do the magic ...

This last line will replay your local commits on top of the fetched
changes and will automatically stash your local changes before pulling
and reapply them after completion to avoid conflicts with uncommitted
changes.

3.2. Merging the updated code into your branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have synched your fork, the ``dev`` branch on your fork matches
the upstream version. You should now get the code on your computer with
the following commands:

::

   cd ramses          # go to your local ramses directory
   git switch dev     # switch to the dev branch
   git pull           # download the updated code

Once you have the updated code, you need to incoroporate it into your
``my-new-feature`` branch. This operation is called a **merge** in Git,
and is done as follows.

::

   git switch my-new-feature
   git merge dev

This will likely generate **conflicts**, which are occurences where Git
cannot decide automatically how combine code from two versions. In such
cases, Git will return ``CONFLICT`` warnings with the list of files that
are concerned. You need to edit the files to resolve the conflicts.
Conflicts are highlighted with the following format:

::

   <<<<<<< HEAD
     ! Here some code in my-new-feature.
     real(dp) :: my_new_feature_secret_variable
   =======
     ! Here is some different code from dev, at the same location...
     real(dp),paramter :: thirteen = 13.0d0
   >>>>>>> dev

In this simple example, you want to keep the two codes and simply have
to remove the lines introduced by Git to produce:

::

     ! Here some code in my-new-feature.
     real(dp) :: my_new_feature_secret_variable

     ! Here is some different code from dev, at the same location...
     real(dp),paramter :: thirteen = 13.0d0

Once you have solved all conflicts, simply commit with

::

   git commit  # Git will suggest an adequate commit message in this case
   git push    # broadcast your commit

.. container:: info

   **Exercise:** Set your upstream to Mentuhotep-II/ramses. Make some
   minor change in the file amr/constants.f90. Wait for a greenlight
   from the lecturer (who will introduce some change in the upstream
   version) and sync your fork with the upstream version.

4. Running the test suite
-------------------------

All changes should be tested to make sure they didn’t break anything.
This can be done by running the test suite as explained in the
`documentation <https://ramses-organisation.readthedocs.io/en/dev/wiki/Testing.html>`__.
Note that they are also run on the github repository before any merge.

.. container:: info

   **Exercise:** Run the all the hydro tests. Make a change and see if
   the test suite still pass. Then explore the ``tests`` folder to
   understand how the tests are done.

5. Creating a pull request
--------------------------

Creating a pull request is the way to incorporate your developement in
the upstream version of the code. The simplest way to do this, after you
committed and pushed your work, is from the GitHub web interface. Simply
click on the “Compare and Create a Pull Request” button and follow the
instructions.

.. raw:: html

   <!-- </div> <!-- perhaps -->
