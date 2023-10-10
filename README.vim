Useful commands:
:set number - see line numbers
:set list - see tabs and white spaces (set nolist to not see special char).
:verbose set ts? et? - to see tabstop and expandtab setting

To replace tabs with spaces:
First, decide how many spaces you want your tabs to be 
converted to. Lets say you want each tab to be 2 spaces. You then do:

:set ts=1  (In my case, I used tab=1 space)
:set et
:%retab!


