use clap::Parser;
use virust_locator::{config::Args, locator};

fn main() {
    let args = Args::parse().validate().unwrap_or_else(|err| {
        eprintln!("{} {}", "\x1b[1;91mError:\x1b[0m", err);
        std::process::exit(1);
    });

    let loc: Vec<Option<locator::Locator>> = locator::Locator::build(&args).unwrap_or_else(|err| {
        eprintln!("{} {}", "\x1b[1;91mError:\x1b[0m", err);
        std::process::exit(1);
    });

    print_loc_vec(loc);
}

fn print_loc_vec(loc: Vec<Option<locator::Locator>>) {
    for l in loc {
        if l.is_none() {
            eprintln!("{} {}", "\x1b[1;91mError:\x1b[0m", "Locator not found");
            std::process::exit(1);
        } else {
            println!("{}", l.unwrap());
        }
    }
}
