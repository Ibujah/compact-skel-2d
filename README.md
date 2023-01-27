# One-step compact skeletonization

This software is an implementation of the compact skeletonization method.

One-step compact skeletonization, Durix B., Morin G., Chambon S., Mari J.-L. and Leonard K., Eurographics 2019

## Website:

http://durix.perso.enseeiht.fr/

## Instructions

Build: 

```
cargo build
```

Run: 

```
cargo run
```

Help: 

```
cargo run -- --help
```

## Few examples

Default test: 
```
cargo run -- --imgfile ressources/rat.png --outfile rat_skel.png
```

Small square 8-connected
```
cargo run -- --imgfile ressources/rat_glitch.png --outfile rat_glitch_skel.png
```

Small square:
```
cargo run -- --imgfile ressources/square_test.png --outfile square_test_skel.png
```

Several connected components:
```
cargo run -- --imgfile ressources/cc_test.png --outfile cc_test_skel.png
```


## Misc

If you find a bug, or have a question, do not hesitate to ask!
