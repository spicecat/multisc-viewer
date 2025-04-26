<script lang="ts">
  import DataTableSearch from "./DataTableSearch.svelte";

  const { publications }: { publications: Publication[] } = $props();

  const publicationData = $derived(
    publications.map((publication) => ({
      ...publication,
      publication: {
        name: publication.name,
        href: `/publication/${publication.publicationId}`,
      },
      year: [...new Set(publication.datasets.map((ds) => ds.year))],
      species: [...new Set(publication.datasets.map((ds) => ds.species))],
      author: [...new Set(publication.datasets.map((ds) => ds.author))],
      disease: [
        ...new Set(publication.datasets.map((ds) => ds.disease).flat()),
      ],
      cellType: [...new Set(publication.datasets.map((ds) => ds.cellType))],
      datasets: publication.datasets,
    }))
  );

  const publicationColumns: Column[] = [
    { key: "publication", label: "Publication", url: true },
    { key: "description", label: "Description" },
    { key: "year", label: "Year" },
    { key: "species", label: "Species" },
    { key: "author", label: "Authors" },
    { key: "disease", label: "Disease" },
    { key: "cellType", label: "Cell Type" },
    { key: "datasets", label: "Datasets" },
  ];
</script>

<DataTableSearch
  label="Publications"
  data={publicationData}
  columns={publicationColumns}
/>
