<script lang="ts">
	import PlotsForm from '$lib/components/Plots/Form.svelte';
	import type { PageProps } from '../$types';

	let { data }: PageProps = $props();
	let { datasets } = $derived(data);
</script>

<svelte:head>
	<title>{data.meta.title}</title>
	<meta name="description" content={data.meta.description} />
</svelte:head>

<section class="mx-auto">
	<h2 class="text-center h2">Datasets</h2>
	<PlotsForm {datasets} />
</section>
