{ nixpkgs, ... }:
let
  pkgs = import nixpkgs {
    system = "x86_64-linux";
  };
in
  with pkgs;
  stdenv.mkDerivation {
    name = "gcut";
    version = "1.0.0";

  # https://nixos.org/nix/manual/#builtin-filterSource
  src = builtins.filterSource
  (path: type: lib.cleanSourceFilter path type
  && baseNameOf path != "matlab") ./.;

  buildInputs = [
    python38Packages.ortools
  ];
}
