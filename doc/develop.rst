Information for developers
============================

The PetClaw repository is hosted on BitBucket using Mercurial
for version control, at http://bitbucket.org/ketch/petclaw.
Some brief development guidelines:

    * Don't force a push!  There should only be one head on the
      Bitbucket repo.
    * PetClaw will never depend on the tip of the development version
      of another package, except for clawpack4petclaw, which is 
      essentially part of the PetClaw project too.  
    * PetClaw depends on several other packages.  Some known working
      builds are list on the wiki.  If you make a change that requires
      an update to one of PetClaw's dependencies, get approval from
      David Ketcheson before pushing it.  In your commit message,
      include the tag or changeset of the version of the updated
      package that you're using.  Also add this information to the
      known working builds.
