def fix(fname):
    im = Image.open(fname)
    newim = im.convert("RGBA")

    source = im.split()
    R, G, B = 0, 1, 2
    rmask = source[R].point(lambda i: i == 255 and 255)
    gmask = source[G].point(lambda i: i > 0 and 255)
    bmask = source[B].point(lambda i: i > 0 and 255)

    out = Image.new("RGBA", im.size, None)

    newim.paste(out, None, rmask)
    newim.paste(im, None, gmask)
    newim.paste(im, None, bmask)

    x, y = newim.size
    for i in range(x / y):
        box = (i * y, 0, y * (i + 1), y)
        region = newim.crop(box)
        region.save("%s_%s.png" % (i, y))  