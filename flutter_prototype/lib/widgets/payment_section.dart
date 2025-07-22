import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../services/payment_service.dart';

class PaymentSection extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Container(
      width: double.infinity, // Full width automatically
      padding: EdgeInsets.all(30),
      color: Colors.black, // Following ui.mdc - no gradients
      child: Column(
        children: [
          PaymentHeader(),
          SizedBox(height: 30),
          PaymentForm(),
        ],
      ),
    );
  }
}

class PaymentHeader extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Row(
      mainAxisAlignment: MainAxisAlignment.spaceBetween,
      children: [
        Expanded(
          child: Column(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Text(
                'Payment Required',
                style: Theme.of(context).textTheme.headlineSmall?.copyWith(
                  fontSize: 24,
                ),
              ),
              SizedBox(height: 8),
              Text(
                '\$0.25 per analysis',
                style: TextStyle(
                  color: Colors.white.withOpacity(0.7), // Secondary text
                  fontSize: 16,
                ),
              ),
            ],
          ),
        ),
        IconButton(
          icon: Icon(Icons.close, color: Colors.white),
          onPressed: () {
            context.read<PaymentService>().hidePayment();
          },
          padding: EdgeInsets.all(8),
          constraints: BoxConstraints(minWidth: 40, minHeight: 40), // Touch-friendly
        ),
      ],
    );
  }
}

class PaymentForm extends StatefulWidget {
  @override
  _PaymentFormState createState() => _PaymentFormState();
}

class _PaymentFormState extends State<PaymentForm> {
  final _nameController = TextEditingController();
  
  @override
  Widget build(BuildContext context) {
    return Column(
      children: [
        // Card input (Stripe widget would go here)
        Container(
          width: double.infinity,
          padding: EdgeInsets.all(16),
          decoration: BoxDecoration(
            color: Colors.white.withOpacity(0.08), // Interactive element
            borderRadius: BorderRadius.circular(8),
          ),
          child: Text(
            'Card Details (Stripe Element)',
            style: TextStyle(color: Colors.white.withOpacity(0.7)),
          ),
        ),
        
        SizedBox(height: 20),
        
        // Name input
        TextField(
          controller: _nameController,
          decoration: InputDecoration(
            labelText: 'Name',
            hintText: 'Optional',
          ),
          style: TextStyle(color: Colors.white),
        ),
        
        SizedBox(height: 20),
        
        // Submit button
        Consumer<PaymentService>(
          builder: (context, paymentService, child) {
            return SizedBox(
              width: double.infinity,
              child: ElevatedButton(
                onPressed: paymentService.isProcessing 
                  ? null 
                  : () => _handlePaymentSetup(context),
                child: paymentService.isProcessing
                  ? Row(
                      mainAxisAlignment: MainAxisAlignment.center,
                      children: [
                        SizedBox(
                          width: 16,
                          height: 16,
                          child: CircularProgressIndicator(
                            strokeWidth: 2,
                            valueColor: AlwaysStoppedAnimation<Color>(Colors.white),
                          ),
                        ),
                        SizedBox(width: 8),
                        Text('Setting up...'),
                      ],
                    )
                  : Text('Setup'),
              ),
            );
          },
        ),
      ],
    );
  }
  
  Future<void> _handlePaymentSetup(BuildContext context) async {
    final paymentService = context.read<PaymentService>();
    
    try {
      await paymentService.setupPayment(
        name: _nameController.text.isEmpty ? null : _nameController.text,
      );
      
      // Success - payment section will auto-hide via state management
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Payment setup successful!')),
      );
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Setup failed: $error')),
      );
    }
  }
  
  @override
  void dispose() {
    _nameController.dispose();
    super.dispose();
  }
} 