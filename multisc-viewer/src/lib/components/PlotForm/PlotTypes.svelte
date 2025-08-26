<script lang="ts">
	import { page } from '$app/state';
	import { PlotType } from '$lib/types/plot';

	let pt = $derived(page.url.searchParams.getAll('pt'));
</script>

<label class="label flex gap-4">
	<span class="label-text">Plots</span>
	{#each Object.entries(PlotType) as [plotTypeKey, plotType] (`plot-type-${plotType}`)}
		<label class="flex items-center space-x-2">
			<input
				type="checkbox"
				name="pt"
				value={plotType}
				class="checkbox"
				checked={pt.includes(plotType)}
			/>
			<p>{plotTypeKey}</p>
		</label>
	{/each}
</label>
