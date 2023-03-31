# One-step compact skeletonization

This software is an implementation of the compact skeletonization method.

One-step compact skeletonization, Durix B., Morin G., Chambon S., Mari J.-L. and Leonard K., Eurographics 2019

## Website:

http://durix.perso.enseeiht.fr/

## Instructions

Build: 

```
cargo build --release
```

Run: 

```
cargo run --release
```

Help: 

```
cargo run --release -- --help
```

## Few examples

Default test: 
```
cargo run --release -- --imgin ressources/rat.png --imgout rat_skel.png
```

Small square 8-connected
```
cargo run --release -- --imgin ressources/rat_glitch.png --imgout rat_glitch_skel.png
```

Small square:
```
cargo run --release -- --imgin ressources/square_test.png --imgout square_test_skel.png
```

Several connected components:
```
cargo run --release -- --imgin ressources/cc_test.png --imgout cc_test_skel.png
```


## Misc

If you find a bug, or have a question, do not hesitate to ask!
