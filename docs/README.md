Welcome to pyEPR documentation:beers:!
===================
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.png?v=103)](https://github.com/zlatko-minev/pyEPR)
[![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/zlatko-minev/pyEPR)
[![star this repo](http://githubbadges.com/star.svg?user=zlatko-minev&repo=pyEPR&style=flat)](https://github.com/zlatko-minev/pyEPR)
[![fork this repo](http://githubbadges.com/fork.svg?user=zlatko-minev&repo=pyEPR&style=flat)](https://github.com/zlatko-minev/pyEPR/fork)


Find the docs at [Read the Docs](http://pyepr-docs.readthedocs.io)

## BUILDING THE DOCS


This folder contains the src and resulting doc files. We use sphinx. 

#### Setup 
Install the theme 

```sh
  pip install sphinx-rtd-theme
```



# First time docs developers

### First time setting up Sphinx 

Notes for developers.

1. Check out this [tut](http://www.ittc.ku.edu/kusp/new/howto/sphinx/setup.html) for quicksetup of sphinx  
``` 
   conda install sphinx numpydoc
   (or pip install -U Sphinx)
```

2. Install read [the docs theme](https://github.com/rtfd/sphinx_rtd_theme) and set in the config `html_theme = "sphinx_rtd_theme"`
``` 
  pip install sphinx_rtd_theme
```

3. Set up autodoc [tut](https://medium.com/@eikonomega/getting-started-with-sphinx-autodoc-part-1-2cebbbca5365) and create full module [doc](https://docs-python2readthedocs.readthedocs.io/en/master/code-doc.html). In the `docs` directory 
```
  sphinx-apidoc -f -o source/ ../pyEPR -o source/api --no-toc -M -e
  make html
```
You can also use this to update the doc tree.

# Updating `readthedocs.org`

The docs get update automatically each time you make a commit to the master branch. To check on the status of the docs build or failing, navigate to `https://readthedocs.org/projects/pyepr-docs/builds/` and select the latest build to see its progress and report. 

To manually initiate a rebuild:  Navigate to `https://readthedocs.org/projects/pyepr-docs/`. Here, if you have write access you can initiate a build. 
