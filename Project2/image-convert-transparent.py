from PIL import Image


def transPNG(srcImageName, dstImageName):
    img = Image.open(srcImageName)
    img = img.convert("RGBA")
    datas = img.getdata()
    newData = list()
    for item in datas:
        if item[0] > 230 and item[1] > 230 and item[2] > 230:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)

    img.putdata(newData)
    img.save(dstImageName, "PNG")


if __name__ == '__main__':
    transPNG("C:\Liang\Googlebox\Conference\ESEB201808\Powertable_all.jpg", "C:\Liang\Googlebox\Conference\ESEB201808\Powertable_all-con.jpg")