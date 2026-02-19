<script lang="ts">
import { ArrowRight } from "@lucide/svelte";
import { Progress } from "@skeletonlabs/skeleton-svelte";
import { navigating, page } from "$app/state";
import DatasetsTable from "$lib/components/Data/DatasetsTable.svelte";
import GenesTable from "$lib/components/Data/GenesTable.svelte";
import type { DatasetsMap, GenesRowsMap } from "$lib/types/daemon";

import Groupings from "./Groupings.svelte";
import PlotTypes from "./PlotTypes.svelte";

const { datasets, genesRows }: { datasets: DatasetsMap; genesRows?: GenesRowsMap } =
	$props();
</script>

<form action="/plots" class="space-y-4">
	<input type="hidden" name="pub" value={page.params.pubId} />
	{#if genesRows}
		<div class="mx-auto flex gap-4">
			<div>
				<GenesTable {datasets} {genesRows} />
			</div>
			<div class="space-y-4">
				<Groupings />
				<PlotTypes />
			</div>
		</div>
	{/if}
	<DatasetsTable {datasets}>
		<button type="submit" class="my-auto btn flex items-center preset-filled-primary-500">
			Plot Datasets
			{#if navigating.to}
				<Progress value={null}>
					<Progress.Circle style="--size: 24px; --thickness: 4px;">
						<Progress.CircleTrack />
						<Progress.CircleRange />
					</Progress.Circle>
				</Progress>
			{:else}
				<ArrowRight size="24" />
			{/if}
		</button>
	</DatasetsTable>
</form>
