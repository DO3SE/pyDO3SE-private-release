import click


@click.group()
def cli():
    """Main cli entrypoint."""
    click.echo("Welcome to do3se_met")


if __name__ == "__main__":
    cli()


def test_cli():
    cli()
