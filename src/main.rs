use mimalloc::MiMalloc;
use clap::Parser;
use clap::Subcommand;

static GLOBAL: MiMalloc = MiMalloc;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
struct Args {
    #[clap(short, long)]
    debug: bool,

    #[clap(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Preprocess {
        list: String,
    }
}

fn main() {
    let _ = Args::parse();
    println!("Hello, world!");
}
