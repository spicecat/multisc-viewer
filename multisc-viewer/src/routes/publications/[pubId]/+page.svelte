<script lang="ts">
	import PlotsForm from '$lib/components/Plots/Form.svelte';
	import Publication from '$lib/components/Data/Publication.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets, publication } = $derived(data);
</script>

<svelte:head>
	<title>{publication.name ?? publication._id}</title>
</svelte:head>

<Publication {publication} />

<hr class="hr" />

<section class="mx-auto">
	<h2 class="text-center h2">Datasets</h2>
	<PlotsForm {datasets} />
</section>
