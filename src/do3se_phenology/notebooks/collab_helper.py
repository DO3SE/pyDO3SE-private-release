# @title Download libraries from github
# @markdown Code Hidden

def setup_dependencies(

):
    try:
        import do3se_phenology
    except ImportError:
        if 'google.colab' in str(get_ipython()):
            print('Running on CoLab')
            print("""
            This notebook requires access to private repositories. To access follow the steps below:
            1. Get an access key from your github account: https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token
            2. Save the key to a file in the following location on your google drive: `My Drive/access/collabaccess.txt`. The file should include your user on the first line and your token on the second.
            """)
            from google.colab import drive
            drive.mount('/content/drive')

            !mkdir - p ~ / .access
            !cp "/content/drive/My Drive/access/collabaccess.txt" ~ / .access / config
            import os
            creds = open(f'{os.path.expanduser("~")}/.access/config')
            creds_parsed = creds.read().splitlines()
            user, token = creds_parsed
            !pip install git + https: // {user}: {token} @ github.com / SEI - DO3SE / thermal_time.git
            !pip install git + https: // {user}: {token} @ github.com / SEI - DO3SE / do3se_phenology.git
            !pip install git + https: // {user}: {token} @ github.com / SEI - DO3SE / do3se_met.git
            creds.close()
            creds = None
            creds_parsed = None
