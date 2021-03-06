[![Build Status](https://travis-ci.org/yageek/lambertgo.svg?branch=master)](https://travis-ci.org/yageek/lambertgo) [![GoDoc](https://godoc.org/github.com/yageek/lambertgo?status.png)](https://godoc.org/github.com/yageek/lambertgo)


#Description
A simple Go package to convert coordinates in Lambert projections system to and from GPS WGS84 coordinates. It is based on the [IGN algorithms and methods](http://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NTG_71.pdf)

# Install
	go get github.com/YaGeek/lambertgo

# Usage

```go
import lgo "github.com/yageek/lambertgo"
import "fmt"

//Declares origin point and translated point
var point *lgo.Point= &lgo.Point{994272.661,113467.422,0,lgo.Meter}

//Converts pointOrg in Lambert Zone 1 to WGS84 - Results in Radian
point.ToWGS84(lgo.LambertI)

//Converts to Degree
point.ToDegree()

fmt.Printf("Latitude:%.5f - Longitude:%.5f",point.Y,point.X)

point.ToRadian()

point.ToLambert(lambertII_e)
```

#License
Copyright (c) 2013 Yannick Heinrich - Released under the GPLv2 License.

