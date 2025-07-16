#!/bin/bash

# Load configuration
source "$(dirname "$0")/config.sh"

echo "🌐 Domain Setup Status for $DOMAIN_NAME"
echo "======================================"
echo ""

# Check DNS Zone Status
echo "📍 DNS Zone Status:"
if gcloud dns managed-zones describe $DNS_ZONE_NAME --format="value(name,dnsName,description)" 2>/dev/null; then
    echo "✅ DNS zone exists"
else
    echo "❌ DNS zone not found"
fi
echo ""

# Check DNS Records
echo "📋 DNS Records:"
if gcloud dns record-sets list --zone=$DNS_ZONE_NAME --format="table(name,type,ttl,rrdatas)" 2>/dev/null; then
    echo "✅ DNS records configured"
else
    echo "❌ Cannot list DNS records"
fi
echo ""

# Check Domain Verification Status
echo "🔍 Domain Verification Status:"
VERIFIED_DOMAIN=$(gcloud domains list-user-verified --filter="domain:$DOMAIN_NAME" --format="value(domain)" 2>/dev/null)
if [ -n "$VERIFIED_DOMAIN" ]; then
    echo "✅ Domain verified: $VERIFIED_DOMAIN"
else
    echo "⏳ Domain verification pending"
    echo "   → Go to: https://search.google.com/search-console/welcome"
    echo "   → Add $DOMAIN_NAME and verify ownership"
fi
echo ""

# Check Domain Mappings
echo "🔗 Domain Mappings:"
MAPPINGS=$(gcloud beta run domain-mappings list --region=$REGION --format="table(name,domain,status)" 2>/dev/null)
if [ -n "$MAPPINGS" ]; then
    echo "$MAPPINGS"
else
    echo "⏳ No domain mappings found"
    echo "   → Run after domain verification:"
    echo "   → gcloud beta run domain-mappings create --service=$FUNCTION_NAME --domain=$DOMAIN_NAME --region=$REGION"
    echo "   → gcloud beta run domain-mappings create --service=$FUNCTION_NAME --domain=www.$DOMAIN_NAME --region=$REGION"
fi
echo ""

# Check Cloud Function Status
echo "☁️ Cloud Function Status:"
FUNCTION_STATUS=$(gcloud functions describe $FUNCTION_NAME --region=$REGION --format="value(name,status,httpsTrigger.url)" 2>/dev/null)
if [ -n "$FUNCTION_STATUS" ]; then
    echo "✅ Cloud Function: $FUNCTION_STATUS"
else
    echo "❌ Cloud Function not found"
fi
echo ""

# Check if nameservers are properly configured
echo "🔧 Nameserver Configuration:"
echo "Current Google Cloud nameservers:"
gcloud dns managed-zones describe $DNS_ZONE_NAME --format="value(nameServers)" 2>/dev/null | tr ';' '\n' | sed 's/^/   /'
echo ""

# Check DNS propagation
echo "🌍 DNS Propagation Check:"
if command -v dig &> /dev/null; then
    echo "Checking $DOMAIN_NAME nameservers:"
    dig $DOMAIN_NAME NS +short 2>/dev/null | sed 's/^/   /' || echo "   ❌ Cannot check nameservers"
    
    echo "Checking www.$DOMAIN_NAME:"
    dig www.$DOMAIN_NAME A +short 2>/dev/null | sed 's/^/   /' || echo "   ❌ Cannot resolve www.$DOMAIN_NAME"
else
    echo "⚠️ dig command not available - install bind-utils to check DNS propagation"
fi
echo ""

# Provide actionable next steps
echo "✅ Next Steps:"
echo "1. Update nameservers at Namecheap to Google Cloud nameservers above"
echo "2. Wait 24-48 hours for DNS propagation"
echo "3. Complete domain verification in Google Search Console"
echo "4. Run domain mapping commands shown above"
echo "5. Test: https://$DOMAIN_NAME and https://www.$DOMAIN_NAME"
echo ""

# Check if we can reach the current function URL
echo "🔗 Current Function URL Test:"
FUNCTION_URL=$(gcloud functions describe $FUNCTION_NAME --region=$REGION --format="value(httpsTrigger.url)" 2>/dev/null)
if [ -n "$FUNCTION_URL" ]; then
    echo "Testing: $FUNCTION_URL"
    if curl -s --max-time 5 "$FUNCTION_URL" > /dev/null 2>&1; then
        echo "✅ Function is responding"
    else
        echo "❌ Function not responding (may be cold start)"
    fi
else
    echo "⚠️ Cannot get function URL - checking function status..."
    gcloud functions describe $FUNCTION_NAME --region=$REGION --format="value(status)" 2>/dev/null || echo "❌ Function not found"
fi 