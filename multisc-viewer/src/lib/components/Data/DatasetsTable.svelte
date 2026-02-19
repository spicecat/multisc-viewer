<script lang="ts">
import { Globe } from "@lucide/svelte";
import { uniq } from "lodash-es";
import type { Snippet } from "svelte";
import { resolve } from "$app/paths";
import { page } from "$app/state";
import Author from "$lib/components/Data/Author.svelte";
import DataTable from "$lib/components/List/DataTable.svelte";
import Search from "$lib/components/List/Search.svelte";
import Tabs from "$lib/components/List/Tabs.svelte";
import type { DatasetsMap } from "$lib/types/daemon";
import OntologyTerm from "./OntologyTerm.svelte";

const { datasets, children }: { datasets: DatasetsMap; children?: Snippet } =
	$props();

let selected = $state(page.url.searchParams.getAll("ds"));

const data = $derived(Object.values(datasets));
let items = $state<typeof data>([]);
let tags = $state<string[]>([]);

const searchOptions = {
	keys: [
		"author.name",
		"citation.name",
		"citation.url",
		"publicationDate",
		"description",
		"healthCondition.displayName",
		"healthCondition.identifier",
		"healthCondition.name",
		"name",
		"species.displayName",
		"species.identifier",
		"species.name",
		"cellType.displayName",
		"cellType.identifier",
		"cellType.name",
		"tissue.displayName",
		"tissue.identifier",
		"tissue.name",
		"url",
	],
};
</script>

<div class="flex gap-4">
  <Search name="Datasets" {data} bind:items bind:tags {searchOptions} />
  {@render children?.()}
</div>
<Tabs
  tabs={uniq(
    data.flatMap((ds) => (ds.cellType ? ds.cellType.map((ot) => ot.name) : [])),
  ).sort()}
  data={items}
  filter={(cell) => (ds) =>
    Boolean(ds.cellType?.some((ct) => cell === ct.name))}
>
  {#snippet panels(cellItems)}
    <DataTable data={cellItems} selected={selected.map((s) => datasets[s])}>
      <tr>
        <th>
          <input
            type="checkbox"
            class="checkbox"
            aria-label="Select all datasets"
            checked={cellItems.every((ds) => selected.includes(ds.id))}
            onchange={(e) =>
              (selected = e.currentTarget.checked
                ? cellItems.map((ds) => ds.id)
                : [])}
          />
          <span class="font-bold">{selected.length}</span>
        </th>
        <th>Dataset</th>
        <th>Cell Type</th>
        <th>Tissue</th>
        <th>Health Condition</th>
        <th>Species</th>
        <th>Author</th>
        <th>Date</th>
      </tr>
      {#snippet row(ds)}
        {@const name = ds.displayName ?? ds.name ?? ds.id}
        <td>
          <input
            type="checkbox"
            class="checkbox"
            aria-label="Select dataset {name}"
            name="ds"
            value={ds.id}
            bind:group={selected}
          />
        </td>
        <td>
          <button type="button" class="btn preset-tonal-primary text-wrap">
            <a
              href={resolve("/datasets/[dsId]", { dsId: ds.id })}
              title={`${name} dataset page`}
            >
              <!-- eslint-disable-next-line svelte/no-at-html-tags -->
              {@html name}
            </a>
            {#if ds.url}
              <!-- eslint-disable svelte/no-navigation-without-resolve -->
              <a
                class="anchor"
                href={ds.url}
                target="_blank"
                title={`${name} dataset source`}
                rel="external noopener noreferrer"
              >
                <Globe />
              </a>
            {/if}
          </button>
        </td>
        <td>
          <OntologyTerm ontologyTerm={ds.cellType} bind:tags tag="cellType" />
        </td>
        <td>
          <OntologyTerm ontologyTerm={ds.tissue} bind:tags tag="tissue" />
        </td>
        <td>
          <OntologyTerm
            ontologyTerm={ds.healthCondition}
            bind:tags
            tag="healthCondition"
          />
        </td>
        <td>
          <OntologyTerm ontologyTerm={ds.species} bind:tags tag="species" />
        </td>
        <td><Author author={ds.author} bind:tags /></td>
        <td class="text-nowrap">{ds.publicationDate}</td>
      {/snippet}
    </DataTable>
  {/snippet}
</Tabs>

<!-- Prevent forgetting selected on page change -->
{#each selected as s, i (`ds-hidden-${i}`)}
  <input type="checkbox" class="hidden" value={s} bind:group={selected} />
{/each}
