# galaxy-tools

Some DeVL/GVL tools for Galaxy.

## Standards
Tools should comply with the [IUC standards](http://galaxy-iuc-standards.readthedocs.io/en/latest/)

**Specifically:**

* [ ] Indentation is correct (4 spaces)
* [ ] Tool version/build ok
* [ ] `<command/>`
  - [ ] Text parameters, input and output files `'single quoted'`
  - [ ] Use of `<![CDATA[ ... ]]>` tags
  - [ ] Parameters of type `text` or having `optional="true"` attribute are checked with `if str($param)` before being used
* [ ] Data parameters have a `format` attribute containing datatypes recognised by Galaxy
* [ ] Tests
  - [ ] Parameters are reasonably covered
  - [ ] Test files are appropriate
* [ ] Help
  - [ ] Valid restructuredText and uses `<![CDATA[ ... ]]>` tags
