package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
	"time"
)

func main() {
	a := os.Args[1]
	if a == "new" {
		newr()
		return
	}
	f, e := ioutil.ReadFile(a)
	fatal(e)
	run(bytes.NewReader(f))
}
func run(r io.Reader) {
	a := Parse(r)
	fmt.Println(a)
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
	E          bool
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
	if g := G[0]; g == 'a' && len(v) == 7 && strings.HasSuffix(v[6], "h") {
		S = flt(strings.TrimSuffix(v[6], "h"))
	} else if g == 'a' && len(v) != 7 {
		pan("a!")
	} else if g == 'b' && len(v) == 7 && v[6] == "e" {
		E = true
	} else if len(v) != 6 {
		pan("len!")
	}
	return R{d, G[0], flt(v[3]), flt(v[4]), flt(v[5]), S, E, 0}
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
