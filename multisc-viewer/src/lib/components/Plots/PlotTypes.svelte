<script lang="ts">
	import { page } from '$app/state';
	import type { PlotsParams } from '$lib/types/daemon';

	const ptOptions: Record<string, PlotsParams['pt'][number]> = {
		UMAP: 'umap',
		Violin: 'vln',
		Feature: 'feat'
	};

	let pt = $state(page.url.searchParams.getAll('pt'));
</script>

<label class="label">
	<span class="label-text">Plots</span>
	<div class="flex gap-4">
		{#each Object.entries(ptOptions) as [ptLabel, p] (`checkbox-plot-type-${p}`)}
			<label class="flex items-center space-x-2">
				<input type="checkbox" name="pt" value={p} class="checkbox" checked={pt.includes(p)} />
				<p>{ptLabel}</p>
			</label>
		{/each}
	</div>
</label>
