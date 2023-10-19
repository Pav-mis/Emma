using ArtifactUtils, Artifacts # Artifacts provides the artifact string macro

julia> add_artifact!(
           "Artifacts.toml",
           "Emma_vertebrate_models",
           "https://github.com/Pav-mis/Emma_vertebrate_models/archive/refs/tags/v1.0.tar.gz",
           force=true,
       )