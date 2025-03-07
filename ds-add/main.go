package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"io/fs"
	"os"
	"os/exec"
	"path"
	"reflect"
	"strconv"
	"strings"

	"golang.org/x/text/cases"
	"golang.org/x/text/language"
)

type CSVDataset struct {
	Year                    int
	Author                  string
	PMID                    string
	Species                 string
	Organ                   string
	BrainStructure          string
	BrainRegion             string
	BrainRegionAbbreviation string
	Disease                 []string
	CellType                string
}

type DSMetadata struct {
	Name    string   `json:"name"`
	Year    int      `json:"year"`
	Region  []string `json:"region"`
	PMID    string
	Species string   `json:"species"`
	Author  string   `json:"author"`
	Disease []string `json:"disease"`
	Size    int64    `json:"size"`
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
	dsRoot := path.Dir(metas[0])
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

	if len(os.Args) > 1 {
		datasetsDir := os.Args[1]

		stat, err := os.Stat(datasetsDir)

		if err != nil {
			fmt.Fprintln(os.Stderr, "Could not find dataset path:", err)
			os.Exit(1)
		}

		if !stat.IsDir() {
			fmt.Fprintln(os.Stderr, fmt.Sprintf("Dataset path '%s' is not a directory:", datasetsDir), err)
			os.Exit(1)
		}

		stat, err = os.Stat(fmt.Sprintf("%s/meta.csv", datasetsDir))

		if err != nil {
			fmt.Fprintln(os.Stderr, "Could not find metadata file:", err)
			os.Exit(1)
		}

		if stat.IsDir() {
			fmt.Fprintln(os.Stderr, "Metadata is not a file")
			os.Exit(1)
		}

		csvMetaFile, err := os.OpenFile(fmt.Sprintf("%s/meta.csv", datasetsDir), os.O_RDONLY, 0)

		if err != nil {
			fmt.Fprintln(os.Stderr, "Could not open metadata file:", err)
			os.Exit(1)
		}

		structFields := reflect.VisibleFields(reflect.TypeFor[CSVDataset]())
		numStructFields := len(structFields)

		reader := csv.NewReader(bufio.NewReader(csvMetaFile))

		fields, err := reader.Read()

		if err == io.EOF {
			fmt.Fprintln(os.Stderr, "No data other than header in CSV")
			os.Exit(1)
		} else if err != nil {
			fmt.Fprintln(os.Stderr, "Error reading header from CSV", err)
			os.Exit(1)
		}

		numFields := len(fields)

		if numFields != numStructFields {
			fmt.Fprintln(os.Stderr, "Incorrect number of headers in CSV")
			os.Exit(1)
		}

		for i := 0; i < numFields; i++ {
			name := fields[i]
			found := false

			for j := 0; j < numFields; j++ {
				if structFields[j].Name == name {
					found = true
					break
				}
			}

			if !found {
				fmt.Fprintf(os.Stderr, "Unknown field '%s' in CSV\n", name)
				os.Exit(1)
			}
		}

		datasets := existing

		for {
			record, err := reader.Read()
			dataset := CSVDataset{}

			if err == io.EOF {
				break
			}

			if err != nil {
				fmt.Fprintln(os.Stderr, "Error reading record from CSV", err)
				os.Exit(1)
			}

			for i := 0; i < numFields; i++ {
				field := reflect.ValueOf(&dataset).Elem().FieldByName(fields[i])

				if field.Type() == reflect.TypeFor[string]() {
					field.Set(reflect.ValueOf(record[i]))
				} else if field.Type() == reflect.TypeFor[int]() {
					val, err := strconv.ParseInt(record[i], 10, int(reflect.TypeFor[int]().Size())*8)

					if err != nil {
						fmt.Fprintln(os.Stderr, "Error parsing integer from CSV", err)
						os.Exit(1)
					}

					field.Set(reflect.ValueOf(int(val)))
				} else if field.Type() == reflect.TypeFor[int64]() {
					val, err := strconv.ParseInt(record[i], 10, 64)

					if err != nil {
						fmt.Fprintln(os.Stderr, "Error parsing integer from CSV", err)
						os.Exit(1)
					}

					field.Set(reflect.ValueOf(val))
				} else if field.Type() == reflect.TypeFor[[]string]() {
					field.Set(reflect.ValueOf(strings.Split(record[i], ", ")))
				}
			}

			dsName := fmt.Sprintf("%s_%d_%s_%s_%s", dataset.Species, dataset.Year, dataset.Author, dataset.BrainRegionAbbreviation, dataset.CellType)

			fmt.Printf("Enter disease for %s:\n", dsName)

			disease, err := bufio.NewReader(os.Stdin).ReadString('\n')
			disease = strings.Trim(disease, "\r\n")

			if err != nil {
				fmt.Fprintln(os.Stderr, "Error reading disease", err)
				os.Exit(1)
			}

			ds := DSMetadata{
				Name: fmt.Sprintf(
					"%s_%d_%s_%s_%s_%s",
					dataset.Species,
					dataset.Year,
					dataset.Author,
					disease,
					dataset.BrainRegionAbbreviation,
					dataset.CellType,
				),
				Year:    dataset.Year,
				PMID:    dataset.PMID,
				Species: dataset.Species,
				Author:  dataset.Author,
				Disease: dataset.Disease,
			}

			dataFileStem := fmt.Sprintf(
				"%s/%s_%d_%s_%s_%s_%s",
				datasetsDir,
				ds.Species,
				ds.Year,
				ds.Author,
				disease,
				strings.ToLower(dataset.BrainRegionAbbreviation),
				dataset.CellType,
			)

			dataFile := fmt.Sprintf("%s_FinalObj_ForWeb.rds", dataFileStem)
			clusterColorFile := fmt.Sprintf("%s_cluster_color.scheme.rds", dataFileStem)
			genotypeColorFile := fmt.Sprintf("%s_genotype_color.scheme.rds", dataFileStem)

			stats, err := os.Stat(dataFile)

			if err != nil {
				fmt.Fprintf(os.Stderr, "Could not find data file for dataset '%s', skipping\n", ds.Name)
			} else {
				ds.Size = stats.Size()

				err = os.Mkdir(fmt.Sprintf("datasets/%s", ds.Name), 0750)

				if err != nil {
					fmt.Fprintln(os.Stderr, "Failed to create dataset dir:", err)
					os.Exit(1)
				}

				os.Rename(dataFile, fmt.Sprintf("datasets/%s/data.rds", ds.Name))
				os.Rename(clusterColorFile, fmt.Sprintf("datasets/%s/cluster.colors.rds", ds.Name))
				os.Rename(genotypeColorFile, fmt.Sprintf("datasets/%s/genotype.colors.rds", ds.Name))

				cwd, err := os.Getwd()

				if err != nil {
					fmt.Fprintln(os.Stderr, "Failed to get CWD:", err)
					os.Exit(1)
				}

				err = os.Chdir(fmt.Sprintf("datasets/%s", ds.Name))

				if err != nil {
					fmt.Fprintln(os.Stderr, "Failed to change CWD:", err)
					os.Exit(1)
				}

				cmd := exec.Command("Rscript", "../../genes.r")
				err = cmd.Start()

				if err != nil {
					fmt.Fprintln(os.Stderr, "Failed to generate genes:", err)
					os.Exit(1)
				}

				err = cmd.Wait()

				if err != nil {
					fmt.Fprintln(os.Stderr, "Failed to generate genes:", err)
					os.Exit(1)
				}

				err = os.Chdir(cwd)

				if err != nil {
					fmt.Fprintln(os.Stderr, "Failed to change CWD:", err)
					os.Exit(1)
				}

				datasets = append(datasets, ds)
			}
		}

		serialized, err := json.MarshalIndent(datasets, "", "\t")

		if err != nil {
			fmt.Fprintln(os.Stderr, "Error serializing to json:", err)
			os.Exit(1)
		}

		meta.WriteAt(serialized, 0)

		fmt.Println("Dataset list updated!")
	} else {
		fmt.Println("Enter new dataset data:")

		ds := DSMetadata{}

		readProp("name", &ds.Name)
		readProp("year", &ds.Year)
		readPropArr("region", &ds.Region)
		readProp("PMID", &ds.PMID)
		readProp("species", &ds.Species)
		readProp("author", &ds.Author)
		readPropArr("disease", &ds.Disease)

		dsPath := path.Join(dsRoot, ds.Name)
		_, err = os.Stat(dsPath)

		if err != nil {
			fmt.Fprintln(os.Stderr, "Could not find dataset path:", err)
			os.Exit(1)
		}

		stat, err := os.Stat(path.Join(dsPath, "data.rds"))

		if err != nil {
			fmt.Fprintln(os.Stderr, "Could not find dataset file:", err)
			os.Exit(1)
		}

		ds.Size = stat.Size()

		_, err = os.Stat(path.Join(dsPath, "genes.json"))

		if err != nil {
			cwd, err := os.Getwd()

			if err != nil {
				fmt.Fprintln(os.Stderr, "Failed to get CWD:", err)
				os.Exit(1)
			}

			err = os.Chdir(dsPath)

			if err != nil {
				fmt.Fprintln(os.Stderr, "Failed to change CWD:", err)
				os.Exit(1)
			}

			cmd := exec.Command("Rscript", "../../genes.r")
			err = cmd.Wait()

			if err != nil {
				fmt.Fprintln(os.Stderr, "Failed to generate genes:", err)
				os.Exit(1)
			}

			err = os.Chdir(cwd)

			if err != nil {
				fmt.Fprintln(os.Stderr, "Failed to change CWD:", err)
				os.Exit(1)
			}
		}

		datasets := append(existing, ds)

		serialized, err := json.MarshalIndent(datasets, "", "\t")

		if err != nil {
			fmt.Fprintln(os.Stderr, "Error serializing to json:", err)
			os.Exit(1)
		}

		meta.WriteAt(serialized, 0)

		fmt.Println("Dataset list updated!")
	}
}

func readProp[T any](name string, loc *T) {
	fmt.Printf("Enter %s: ", name)

	if reflect.TypeOf(loc) == reflect.TypeFor[*string]() {
		rd := bufio.NewReader(os.Stdin)
		str, err := rd.ReadString('\n')

		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading %s: %s\n", name, err.Error())
			os.Exit(1)
		}

		*loc = interface{}(strings.TrimSpace(str)).(T)
	} else {
		_, err := fmt.Scan(loc)

		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading %s: %s\n", name, err.Error())
			os.Exit(1)
		}
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
