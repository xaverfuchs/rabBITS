library(devtools)

create_package("/Users/b1082752/Nextcloud/GitHub repositories/rabBITS") #done

use_github()

use_r()

load_all()


usethis::use_pkgdown()


pkgdown::build_site()



#sticker creation
library(hexSticker)
imgurl <- "logo/rabbit.png"
sticker(imgurl, package="rabBITS", s_x=1, s_y=1, s_width=1,
        p_size=33, p_y = 0.75, p_color="magenta", p_family="sans", p_fontface="bold",
        url="xaverfuchs.github.io/raBITS", u_size = 5, u_color = "white",
        h_color="green",
        white_around_sticker = T,
        filename="logo/logoimg.png")



usethis::use_logo(img = "logo/logoimg.png")



usethis::use_github()

usethis::use_pkgdown_github_pages()



