use clap::Parser;
use virust_locator::{config::Args, locator};

fn main() {
    let args = Args::parse().validate().unwrap_or_else(|err| {
        eprintln!("{} {}", "\x1b[1;91mError:\x1b[0m", err);
        std::process::exit(1);
    });

    let loc = locator::Locator::build(&args).unwrap_or_else(|err| {
        eprintln!("{} {}", "\x1b[1;91mError:\x1b[0m", err);
        std::process::exit(1);
    });

    if loc.is_none() {
        eprintln!("{} {}", "\x1b[1;91mError:\x1b[0m", "Locator not found");
        std::process::exit(1);
    } else {
        println!("{}", loc.unwrap());
    }
}
