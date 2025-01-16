package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io/fs"
	"os"
	"reflect"
	"strings"

	"golang.org/x/text/cases"
	"golang.org/x/text/language"
)

type DSMetadata struct {
	Name    string   `json:"name"`
	Year    int      `json:"year"`
	Region  []string `json:"region"`
	PMID    string
	Species string   `json:"species"`
	Author  string   `json:"author"`
	Disease []string `json:"disease"`
}

func main() {
	fmt.Println("Reading existing datasets...")

	fsys := os.DirFS(".")
	metas, _ := fs.Glob(fsys, "**/meta.json")

	if len(metas) == 0 {
		fmt.Println("Could not find meta.json")
		os.Exit(1)
	}

	// TODO: prioritize deeper-nested meta.json, rather than first
	meta, err := os.OpenFile(metas[0], os.O_RDWR, 0)

	if err != nil {
		fmt.Printf("Error opening file %s: %s\n", metas[0], err.Error())
		os.Exit(1)
	}

	rawExisting := make([]byte, 0)

	buf := make([]byte, 1000)
	for {
		n, _ := meta.Read(buf)

		if n == 0 {
			break
		} else {
			rawExisting = append(rawExisting, buf[0:n]...)
		}
	}

	var existing []DSMetadata
	err = json.Unmarshal(rawExisting, &existing)

	if err != nil {
		fmt.Printf("Error reading existing datasets: %s\n", err.Error())
		os.Exit(1)
	}

	fmt.Println("Enter new dataset data:")

	ds := DSMetadata{}

	readProp("name", &ds.Name)
	readProp("year", &ds.Year)
	readPropArr("region", &ds.Region)
	readProp("PMID", &ds.PMID)
	readProp("species", &ds.Species)
	readProp("author", &ds.Author)
	readPropArr("disease", &ds.Disease)

	datasets := append(existing, ds)

	serialized, err := json.MarshalIndent(datasets, "", "\t")

	if err != nil {
		fmt.Fprintln(os.Stderr, "Error serializing to json:", err)
		os.Exit(1)
	}

	meta.WriteAt(serialized, 0)

	fmt.Println("Dataset list updated!")
}

func readProp[T any](name string, loc *T) {
	fmt.Printf("Enter %s: ", name)
	_, err := fmt.Scan(loc)

	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading %s: %s\n", name, err.Error())
		os.Exit(1)
	}
}

func readPropArr[T any](name string, loc *[]T) {
	arr := make([]T, 0)

	fmt.Printf("Enter %ss (blank to end): \n", name)
	i := 1
	for ; ; i++ {
		var in T

		fmt.Printf("%s %d: ", cases.Title(language.English).String(name), i)
		rd := bufio.NewReader(os.Stdin)
		str, err := rd.ReadString('\n')

		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading element of %s: %s\n", name, err.Error())
			os.Exit(1)
		}

		str = strings.TrimSpace(str)
		if reflect.TypeOf(in) == reflect.TypeFor[string]() {
			if str == "" {
				*loc = arr
				return
			} else {
				arr = append(arr, interface{}(str).(T))
			}
		} else {
			_, err := fmt.Sscan(str, &in)

			if err != nil {
				fmt.Fprintf(os.Stderr, "Error reading element of %s: %s\n", name, err.Error())
				os.Exit(1)
			}

			if reflect.ValueOf(in) == reflect.Zero(reflect.TypeOf(in)) {
				*loc = arr
				return
			} else {
				arr = append(arr, in)
			}
		}
	}
}
