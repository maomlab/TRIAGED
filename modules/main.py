from modules.misc_utils.graphics_utilis import make_banner

def main():
    banner, banner_foot = make_banner()
    print(banner + banner_foot + "\033[0m")
    

if __name__ == "__main__":
    main()