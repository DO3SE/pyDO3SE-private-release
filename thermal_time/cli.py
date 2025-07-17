import click


@click.group()
def cli():
    """Main cli entrypoint."""
    click.echo("Welcome to do3se_phenology")

@cli.command()
def convert_phenology():
    raise NotImplementedError()

if __name__ == "__main__":
    cli()