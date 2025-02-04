SPARTA MCFX-6DoFconv
---

This plugin is a 6DoF convolver based on SPARTA's [6DoFConv](https://github.com/leomccormack/SPARTA?tab=readme-ov-file#plug-in-descriptions), which replaces its convolution engine with that of the [MCFX Convolver](https://github.com/kronihias/mcfx) plugin.

MCFX's convolver uses non-uniform partitioned convolutions that proved much more efficient than the original implementation. 
This enabled to:
1. Reduce audio latency introduced by SPARTA 6DoFConv on buffer sizes <512 (reduced to 0 samples).
2. Reduce greatly CPU usage, reflected on greatly increased real-time ratio.
3. Reduced IR matrix change delay from 2 buffer-sizes to 0.
4. Enabled automatic IR resampling on load (with warning message).

Plus, the source code is now ready for future addition of the following features:
- Multisource Support (not yet implemented)
- Longer crossfade time (not yet implemented)
