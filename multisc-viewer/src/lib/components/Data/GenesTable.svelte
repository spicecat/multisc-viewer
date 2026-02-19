<script lang="ts">
import { Dna } from "@lucide/svelte";
import { page } from "$app/state";
import DataTable from "$lib/components/List/DataTable.svelte";
import Search from "$lib/components/List/Search.svelte";
import type { DatasetsMap, GenesRowsMap } from "$lib/types/daemon";

const { datasets, genesRows }: { datasets: DatasetsMap; genesRows: GenesRowsMap } =
	$props();

let selected = $state(page.url.searchParams.getAll("gene"));

const data = $derived(Object.values(genesRows));
let items = $state<typeof data>([]);
let tags = $state<string[]>([]);

const searchOptions = { 
	keys: [
		"id", 
		"ds", 
		"degs",
	], 
};
</script>

<Search
	name="Differentially Expressed Genes"
	{data}
	bind:items
	bind:tags
	{searchOptions}
/>
<DataTable data={items} selected={data.filter((g) => selected.includes(g.id))}>
	<tr>
		<th>
			<label class="flex items-center space-x-2">
				<input
					type="checkbox"
					class="checkbox"
					aria-label="Select all genes"
					checked={items.every((gene) => selected.includes(gene.id))}
					onchange={(e) => (selected = e.currentTarget.checked ? items.map((gene) => gene.id) : [])}
				/>
				<p class="font-bold">{selected.length}</p>
			</label>
		</th>
		<th>Gene</th>
		{#each Object.keys(datasets) as ds (`genes-th-ds-${ds}`)}
			{@const degs = datasets[ds].deg}
			{#if degs}
				{#each Object.values(degs) as deg, i (`degs-th-deg-${deg.id}`)}
					<th class={i === 0 ? 'border-l border-gray-300' : undefined}>
						{deg.name ?? deg.id}
						<p class="font-bold">
							<span>{datasets[ds].genes} genes</span>
							<span class="text-red-500">{deg.genes} DEGs</span>
						</p>
					</th>
				{/each}
			{:else}
				<th class="border-x border-gray-300">
					{datasets[ds]?.displayName ?? ds}
					<p class="font-bold">
						<span>{datasets[ds].genes} genes</span>
					</p>
				</th>
			{/if}
		{/each}
		<th># of Datasets</th>
	</tr>
	{#snippet row(gene)}
		<td>
			<input
				type="checkbox"
				class="checkbox"
				aria-label="Select gene {gene.id}"
				name="gene"
				value={gene.id}
				bind:group={selected}
			/>
		</td>
		<td>{gene.id}</td>
		{#each Object.keys(datasets) as ds (`genes-td-ds-${ds}`)}
			{@const degs = datasets[ds].deg}
			{#if degs}
				{#each Object.values(degs) as deg, i (`degs-td-deg-${deg.id}`)}
					<td class={i === 0 ? 'border-l border-gray-300' : undefined}>
						<div class="flex justify-center">
							{#if gene.deg.includes(deg.id)}
								<Dna color="red" />
							{:else if gene.gene.includes(ds)}
								<Dna />
							{/if}
						</div>
					</td>
				{/each}
			{:else}
				<td class="border-x border-gray-300">
					<div class="flex justify-center">
						{#if gene.gene.includes(ds)}
							<Dna />
						{/if}
					</div>
				</td>
			{/if}
		{/each}
		<td class="fold-bold text-center">
			<p class="font-bold">
				<span>{gene.gene.length}, </span>
				<span class="text-red-500">{gene.deg.length}</span>
			</p>
		</td>
	{/snippet}
</DataTable>

<!-- Prevent resetting selected on page change -->
{#each selected as s, i (`gene-hidden-${i}`)}
	<input type="checkbox" class="hidden" value={s} bind:group={selected} />
{/each}
