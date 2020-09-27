package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"sort"
	"strconv"
	"strings"
	"time"
)

func main() {
	var g string
	var nw bool
	var s, a bool
	flag.BoolVar(&nw, "new", false, "new")
	flag.StringVar(&g, "g", "abcn", "filter group")
	flag.BoolVar(&s, "s", false, "stats")
	flag.BoolVar(&a, "a", false, "average")
	flag.Parse()
	if nw {
		newr()
		return
	}
	f, e := ioutil.ReadFile("c:/k/p")
	fatal(e)
	l := filter(g, Parse(bytes.NewReader(f)))
	if a {
		l = average(l)
	}
	if s {
		stats(l)
	} else {
		fmt.Println(l)
	}
}
func fatal(e error) {
	if e != nil {
		panic(e)
	}
}

type R struct {
	Date       time.Time
	G          byte
	H, L, B, S float64
	E, Z       bool
	Rep        int
}
type L []R

func (a L) String() string {
	v := make([]string, len(a))
	for i, r := range a {
		v[i] = r.String()
	}
	return strings.Join(v, "\n") + "\n"
}
func (r R) String() string {
	s := "   "
	if r.G == 'a' {
		s = fmt.Sprintf("%2dh", int(r.S))
		if r.Z {
			s = s[:len(s)-1] + "z"
		}
	}
	return fmt.Sprintf("%s %s %c %03d %02d %02d %s %d", r.Date.Weekday().String()[:3], r.Date.Format("2006.01.02"), r.G, int(r.H), int(r.L), int(r.B), s, r.Rep)
}
func Parse(r io.Reader) (a L) {
	sc := bufio.NewScanner(r)
	for sc.Scan() {
		t := sc.Text()
		if len(t) == 0 {
			continue
		}
		v := strings.Fields(sc.Text())
		if n := len(v); n == 0 || n == 3 {
			continue
		} else if n == 6 || n == 7 {
			a = append(a, parse(t))
		} else {
			panic(fmt.Errorf("parse: %q", t))
		}
	}
	fatal(sc.Err())
	return reps(a)
}
func parse(s string) R {
	pan := func(a string) { fatal(fmt.Errorf("%s %s", a, s)) }
	v := strings.Fields(s)
	d, e := time.Parse("Mon 2006.01.02", v[0]+" "+v[1])
	fatal(e)
	if w := d.Weekday().String(); w[:3] != v[0] {
		pan("weekday!")
	}
	if len(v[2]) != 1 || strings.Index("abcn", v[2]) == -1 {
		pan("group!")
	}
	G := []byte(v[2])
	S := 0.0
	E := false
	Z := false
	if g := G[0]; g == 'a' && len(v) == 7 {
		if s := v[6]; strings.HasSuffix(s, "h") || strings.HasSuffix(s, "z") {
			S = flt(s[:len(s)-1])
			if strings.HasSuffix(s, "z") {
				Z = true
			}
		} else {
			panic("a!")
		}
	} else if g == 'a' && len(v) != 7 {
		pan("a!")
	} else if g == 'b' && len(v) == 7 && v[6] == "e" {
		E = true
	} else if len(v) != 6 {
		pan("len!")
	}
	return R{d, G[0], flt(v[3]), flt(v[4]), flt(v[5]), S, E, Z, 0}
}
func flt(s string) float64 {
	f, e := strconv.ParseFloat(s, 64)
	fatal(e)
	return f
}
func reps(a L) L {
	for i, r := range a[1:] {
		if r.Date == a[i].Date && r.G == a[i].G {
			a[1+i].Rep = a[i].Rep + 1
		}
	}
	return a
}
func newr() {
	t := time.Now()
	o3 := func() {
		s := t.Format("Mon 2006.01.02")
		fmt.Printf("%s a\n%s b\n%s c\n\n", s, s, s)
	}
	o3()
	t = t.AddDate(0, 0, 1)
	o3()
	t = t.AddDate(0, 0, 1)
	o3()
}
func filter(g string, a L) (r L) {
	for _, x := range a {
		if strings.IndexByte(g, x.G) != -1 {
			r = append(r, x)
		}
	}
	return r
}
func average(a L) (r L) {
	var x R
	for i, y := range a {
		if y.Rep == 0 {
			if i != 0 {
				r = append(r, x)
			}
			x = y
		} else if y.Rep == 1 {
			x = merge(x, y)
		} else {
			panic("rep")
		}
	}
	return append(r, x)
}
func merge(x, y R) R {
	f := func(a, b float64) float64 {
		return float64((int(a) + int(b)) / 2)
	}
	x.H = f(x.H, y.H)
	x.L = f(x.L, y.L)
	x.B = f(x.B, y.B)
	return x
}
func floats(a L, f func(x R) float64) []float64 {
	var v []float64
	for _, r := range a {
		if x := f(r); x != 0 {
			v = append(v, x)
		}
	}
	return v
}
func st(name string, v []float64) {
	N := len(v)
	if N == 0 {
		return
	}
	Min, Max := v[0], v[0]
	for _, x := range v {
		if x < Min {
			Min = x
		}
		if x > Max {
			Max = x
		}
	}
	b := make([]float64, len(v))
	copy(b, v)
	sort.Float64s(b)
	Med := b[N/2]
	fmt.Printf("%s %3v [%3v %3v] #%d\n", name, Med, Min, Max, N)
}
func stats(a L) {
	st("H", floats(a, func(x R) float64 { return x.H }))
	st("L", floats(a, func(x R) float64 { return x.L }))
	st("B", floats(a, func(x R) float64 { return x.B }))
	st("S", floats(a, func(x R) float64 { return x.S }))
}
