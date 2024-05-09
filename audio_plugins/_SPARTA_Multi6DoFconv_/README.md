SPARTA Multi6DoFconv
---

This plugin is a 6DoF convolver made to extend SPARTA's [6DoFConv](https://github.com/leomccormack/SPARTA?tab=readme-ov-file#plug-in-descriptions) capabilities.

The final goal is to have MIMO convolution (multiple input, multiple output) while 6DofConv is SIMO (only mono input, and multiple output).

The process involved the substitution of SPARTA's convolver engine `tvconv` with the [MCFX Convolver](https://github.com/kronihias/mcfx) engine.

The process brought a few additional features.

The current features ADDED and MISSING are the following:
```diff
+ MCFX Non-uniform partitioned Convolver engine
+ Automatic IR resampling on load (with warning message)
- Multisource Support (not yet implemented)
- Longer crossfade time (not yet implemented)
```