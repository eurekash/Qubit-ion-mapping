{ sources ? import ./nix/sources.nix { }
, pkgs ? import sources.nixpkgs { }
, gitignore ? (import sources.gitignore { inherit (pkgs) lib; }).gitignoreSource
}:

pkgs.stdenv.mkDerivation {
    pname = "qubit-ion-mapping";
    version = "unstable";

    src = gitignore ./.;

    buildPhase = ''
        g++ -o optimizer ./optimizer.cpp -O3
    '';

    installPhase = ''
        mkdir -p $out/bin
        cp optimizer $out/bin/
    '';

    checkPhase = ''
        $out/bin/optimizer ./sample-input.txt 5
    '';
}