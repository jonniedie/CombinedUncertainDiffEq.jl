# CombinedUncertainDiffEq

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jonniedie.github.io/CombinedUncertainDiffEq.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jonniedie.github.io/CombinedUncertainDiffEq.jl/dev)
[![Build Status](https://github.com/jonniedie/CombinedUncertainDiffEq.jl/workflows/CI/badge.svg)](https://github.com/jonniedie/CombinedUncertainDiffEq.jl/actions)
[![Coverage](https://codecov.io/gh/jonniedie/CombinedUncertainDiffEq.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jonniedie/CombinedUncertainDiffEq.jl)

This package attempts to provide a generic utilities for propagating uncertainty through differential equations simulations using DiffEqUncertainty.jl and ReachabilityAnalysis.jl. The main goal is to experiment with different methods for training and validating neural control systems using full envelope uncertainty propagation instead of single-trajectory Monte Carlo approaches.

*Disclaimer*: Don't expect much. I have no idea what I'm doing at this point. This is just an idea that sounded interesting to me so I'm seeing if it will work.