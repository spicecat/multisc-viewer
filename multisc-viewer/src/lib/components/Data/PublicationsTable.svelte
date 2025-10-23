<script lang="ts">
	import { resolve } from '$app/paths';
	import type { Datasets, Publications } from '$lib/types/daemon';
	import Search from '../List/Search.svelte';
	import Chip from '../List/Chip.svelte';

	let { datasets, publications }: { datasets: Datasets; publications: Publications } = $props();

	let data = $derived(
		Object.values(publications).map((pub) => ({
			...pub,
			datasets: pub.datasets.map((ds) => datasets[ds])
		}))
	);
	let items: typeof data = $state([]);
	let tags: string[] = $state([]);

	const searchOptions = {
		keys: [
			'author.name',
			'date',
			'description',
			'journalName',
			'identifier',
			'name',
			'url',
			'datasets._id',
			'datasets.displayName',
			'datasets.name',
			'datasets.healthCondition.displayName',
			'datasets.healthCondition.identifier',
			'datasets.healthCondition.name',
			'datasets.species.displayName',
			'datasets.species.identifier',
			'datasets.species.name',
			'datasets.cellType.displayName',
			'datasets.cellType.identifier',
			'datasets.cellType.name',
			'datasets.tissue.displayName',
			'datasets.tissue.identifier',
			'datasets.tissue.name'
		]
	};
</script>

<Search name="Publications" {data} bind:items bind:tags {searchOptions} />
{#each items as pub (`pub-${pub._id}`)}
	{@const name = pub.name ?? pub._id}
	<div>
		<a
			href={resolve('/publications/[pubId]', { pubId: pub._id })}
			class="anchor text-lg font-medium"
			title={`${name} publication page`}
			rel="external noopener noreferrer"
		>
			<!-- eslint-disable-next-line svelte/no-at-html-tags -->
			{@html name}
		</a>
		<div>
			<span class="block">
				{pub.author?.map((a) => a.name).join(', ')}
			</span>
			<span class="inline-block">
				{pub.journalName}
				{pub.date}
				<!-- eslint-disable svelte/no-navigation-without-resolve -->
				{#if pub.url}
					<a
						href={pub.url}
						class="anchor"
						target="_blank"
						title={`${name} publication source`}
						rel="external noopener noreferrer"
					>
						{pub.identifier ?? pub.url}
					</a>
				{/if}
				<!-- eslint-enable svelte/no-navigation-without-resolve -->
			</span>
			<span class="block space-y-1 space-x-1">
				{#each pub.datasets as ds (`ds-${pub._id}-${ds._id}`)}
					{@const name = ds.displayName ?? ds.name ?? ds._id}
					<Chip bind:tags tag={`datasets=${name}`}>{name}</Chip>
				{/each}
			</span>
		</div>
		<div>
			<p class="description">
				<!-- eslint-disable-next-line svelte/no-at-html-tags -->
				{@html pub.description}
			</p>
		</div>
	</div>
	<hr class="hr" />
{/each}

<style>
	.description {
		display: -webkit-box;
		-webkit-box-orient: vertical;
		-webkit-line-clamp: 3;
		line-clamp: 3;
		overflow: hidden;
	}
</style>
