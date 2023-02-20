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
cargo run --release -- --imgfile ressources/rat.png --outfile rat_skel.png
```

Small square 8-connected
```
cargo run --release -- --imgfile ressources/rat_glitch.png --outfile rat_glitch_skel.png
```

Small square:
```
cargo run --release -- --imgfile ressources/square_test.png --outfile square_test_skel.png
```

Several connected components:
```
cargo run --release -- --imgfile ressources/cc_test.png --outfile cc_test_skel.png
```


## Misc

If you find a bug, or have a question, do not hesitate to ask!
