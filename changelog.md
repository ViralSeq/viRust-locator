# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [0.1.4] - 2025-06-11

### Added

- `Clone` and `PartialEq` trait for struct `Locator`

- `prelude` module

## [0.1.3] - 2025-05-08

### Added

- Documentations

## [0.1.1] - 2025-05-02

### Added

- `algorithm1` now accepts `&[u8]` for better performance under Rayon

### Fixed

- threading error by switching to `Box<dyn Error + Send + Sync>`

## [0.1.0] - 2025-04-30

### Added

- initial build
