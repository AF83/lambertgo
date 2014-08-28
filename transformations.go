// Package helping to transform coordinates from Lambert system into the WGS84 system
package lambertgo

import (
	"fmt"
	"math"
)

// latitudeISOFromLatitude implements ALG0001
func latitudeISOFromLatitude(lat float64, e float64) float64 {
	return math.Log(math.Tan(math.Pi/4+lat/2) * math.Pow((1-e*math.Sin(lat))/(1+e*math.Sin(lat)), e/2))
}

// latitudeISOFromLatitude implement ALG0002
func latitudeFromLatitudeISO(lat_iso float64, e float64, eps float64) float64 {

	phi_0 := 2*math.Atan(math.Exp(lat_iso)) - math.Pi/2
	phi_i := 2*math.Atan(math.Pow((1+e*math.Sin(phi_0))/(1-e*math.Sin(phi_0)), e/2)*math.Exp(lat_iso)) - math.Pi/2

	delta := 100.0
	for delta > eps {
		phi_0 = phi_i
		phi_i = 2*math.Atan(math.Pow((1+e*math.Sin(phi_0))/(1-e*math.Sin(phi_0)), e/2.0)*math.Exp(lat_iso)) - math.Pi/2
		delta = math.Abs(phi_i - phi_0)
	}
	return phi_i
}

// lambertToGeographic implements
func (pt *Point) lambertToGeographic(zone Zone, lon_merid float64, e float64, eps float64) {

	n := lambertN[zone]
	C := lambertC[zone]
	x_s := lambertXs[zone]
	y_s := lambertYs[zone]
	x := pt.X
	y := pt.Y

	R := math.Sqrt(((x-x_s)*(x-x_s) + (y-y_s)*(y-y_s)))
	gamma := math.Atan((x - x_s) / (y_s - y))

	lon := lon_merid + gamma/n
	lat_iso := -1 / n * math.Log(math.Abs(R/C))

	lat := latitudeFromLatitudeISO(lat_iso, e, eps)

	pt.X = lon
	pt.Y = lat
}

// lambertNormal implement ALG0021
func lambertNormal(lat float64, a float64, e float64) float64 {

	sina := math.Sin(lat)
	N := a / math.Sqrt(1-e*e*sina*sina)
	return N
}

func (pt *Point) geographicToCartesian(a float64, e float64) {

	lat := pt.Y
	lon := pt.X
	he := pt.Z

	N := lambertNormal(lat, a, e)

	pt.X = (N + he) * math.Cos(lat) * math.Cos(lon)
	pt.Y = (N + he) * math.Cos(lat) * math.Sin(lon)
	pt.Z = (N*(1-e*e) + he) * math.Sin(lat)
}

func (pt *Point) cartesianToGeographic(meridien float64, a float64, e float64, eps float64) {

	x := pt.X
	y := pt.Y
	z := pt.Z

	lon := meridien + math.Atan(y/x)

	module := math.Sqrt(x*x + y*y)

	phi_0 := math.Atan(z / (module * (1 - (a*e*e)/math.Sqrt(x*x+y*y+z*z))))

	phi_i := math.Atan(z / module / (1 - a*e*e*math.Cos(phi_0)/(module*math.Sqrt(1-e*e*math.Sin(phi_0)*math.Sin(phi_0)))))

	delta := 100.0
	for delta > eps {

		phi_0 = phi_i
		phi_i = math.Atan(z / module / (1 - a*e*e*math.Cos(phi_0)/(module*math.Sqrt(1-e*e*math.Sin(phi_0)*math.Sin(phi_0)))))
		delta = math.Abs(phi_i - phi_0)
	}
	he := module/math.Cos(phi_i) - a/math.Sqrt(1-e*e*math.Sin(phi_i)*math.Sin(phi_i))

	pt.X = lon
	pt.Y = phi_i
	pt.Z = he
	pt.Unit = Radian
}

// ToWGS84 converts coordinates expressed in Meter in the lambert system to Radian in the WGS84 system.
// It takes the lambert Zone ine parameters
func (pt *Point) ToWGS84(zone Zone) {

	if pt.Unit != Meter {
		fmt.Errorf("Could not transform Point which is not in METER\n")
		return
	}
	if Lambert93 == zone {
		pt.lambertToGeographic(zone, IERSLongitudeMeridian, EWGS84, DefaultEPS)
		pt.Unit = Radian
	} else {
		pt.lambertToGeographic(zone, ParisLongitudeMeridian, EClarkIGN, DefaultEPS)
		pt.Unit = Radian
		pt.geographicToCartesian(AClarkIGN, EClarkIGN)

		pt.X -= 168
		pt.Y -= 60
		pt.Z += 320

		pt.cartesianToGeographic(GreenwichLongitudeMeridian, AWGS84, EWGS84, DefaultEPS)

	}

}

/*
 *  http://geodesie.ign.fr/contenu/fichiers/documentation/pedagogiques/TransformationsCoordonneesGeodesiques.pdf
 *  3.4 Coordonnées géographiques Lambert
 */
func (pt *Point) geographicToLambert(zone Zone, lon_merid float64) {
	long := pt.X
	lat := pt.Y
	n := lambertN[zone]
	c := lambertC[zone]
	e := EWGS84
	xs := lambertXs[zone]
	ys := lambertYs[zone]

	lat_iso := latitudeISOFromLatitude(lat, e)

	pt.X = xs + c*math.Exp(-n*lat_iso)*math.Sin(n*(long-lon_merid))
	pt.Y = ys - c*math.Exp(-n*lat_iso)*math.Cos(n*(long-lon_merid))
	pt.Unit = Meter

}

// ToLambert implements WGS84 -> Lambert2_E
func (pt *Point) ToLambert(zone Zone) {

	if pt.Unit != Radian {
		fmt.Errorf("Could not transform Point which is not in RADIAN\n")
		return
	}
	if Lambert93 == zone {
		fmt.Errorf("NOPE, Lambert 93 is not supported yet.\n")
	} else {
		pt.X = pt.X - GreenwichLongitudeMeridian
		pt.geographicToCartesian(AWGS84, EWGS84)
		pt.X += 168
		pt.Y += 60
		pt.Z -= 320
		pt.cartesianToGeographic(ParisLongitudeMeridian, AWGS84, EWGS84, DefaultEPS)
	}

	pt.geographicToLambert(zone, ParisLongitudeMeridian)

}

// ToLambertAlg003 implements Alg003, pretty useless
func (pt *Point) ToLambertAlg003(zone Zone) {
	if pt.Unit != Radian {
		fmt.Errorf("Could not transform Point which is not in RADIAN. See Point.ToRadian\n")
	}

	long := pt.X
	lat := pt.Y
	n := lambertN[zone]
	c := lambertC[zone]
	e := EClarkIGN
	lambda_c := 0.04079234433 //2.337229167 * math.Pi / 180 //fixme
	xs := lambertXs[zone]
	ys := lambertYs[zone]

	lat_iso := latitudeISOFromLatitude(lat, e)

	x := xs + c*math.Exp(-n*lat_iso)*math.Sin(n*(long-lambda_c))
	y := ys - c*math.Exp(-n*lat_iso)*math.Cos(n*(long-lambda_c))

	pt.X = x
	pt.Y = y
	pt.Unit = Meter
}
